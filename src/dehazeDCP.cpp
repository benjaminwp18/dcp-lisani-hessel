/**
 *
 * Copyright (c) 2023, Jose-Luis Lisani (minor changes), joseluis.lisani@uib.es
 * Based on implementation by He Zhang (https://github.com/He-Zhang/image_dehaze)
 *
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "library/parser.h"
#include "library/libImageFormats.h"

#include <vector>
#include <algorithm>
using namespace std;

#define clip(x, low, high) (x < low)?(low):((x > high)?(high):(x))


//Create auxiliary image for boundary handling
//Use convention dcb|abcdefgh|gfe
//Input parameters:
//  - input image
//  - width and height of input image
//  - r: size of the boundary around the original image
//Output:  image of size width+2*r x height+2*r,
//         the original image starts at position (r, r) and
//         the pixels at the added boundary are obtained by reflection
//         of the original values using the convention dcb|abcdefgh|gfe
double *boundaryHandlingImage(double *in, int width, int height, int r)
{
    double *out = new double[(width+2*r)*(height+2*r)];
    int widthb = width+2*r;

    //interior of output image (Rectangle from pixel (r, r) to (r+width-1, r+height-1)) = input image
    for (int j=0; j < height; j++)
        for (int i=0; i < width; i++)
            out[(i+r)+(j+r)*widthb] = in[i+j*width];

    //Apply reflection dcb|abcdefgh|gfe

    //left side
    for (int j=0; j < height; j++)
        for (int k=1; k <= r; k++)
            out[(r-k)+(j+r)*widthb] = out[(r+k)+(j+r)*widthb];

    //right side
    for (int j=0; j < height; j++)
        for (int k=1; k <= r; k++)
            out[(r+width-1+k)+(j+r)*widthb] = out[(r+width-1-k)+(j+r)*widthb];

    //top
    for (int i=0; i < widthb; i++)
        for (int k=1; k <= r; k++)
            out[i+(r-k)*widthb] = out[i+(r+k)*widthb];

    //bottom
    for (int i=0; i < widthb; i++)
        for (int k=1; k <= r; k++)
            out[i+(r+height-1+k)*widthb] = out[i+(r+height-1-k)*widthb];


    return out;
}


//This function implements line 1 of Algorithm 5 in accompanying IPOL article
//compute min_{over neighborhood} min_{over channels}
//Input parameters:
//  - input image
//  - width, height and number of channels of input image
//  - sr: radius of the square patches used to compute the dark channel
//  - udcp: boolean flag, if TRUE compute dark channel as minimum over G and B values,
//                        else compute it as minimum over R, G and B
//Output:  dark channel image
double *getDarkChannel(double *input, int width, int height, int nchannels, int sr, bool udcp)
{
    int N = width*height;
    double *Jdark = new double[N];

    //compute image of minimum values over R, G and B
    double *minRGB = new double[N];
    if (nchannels == 1) memcpy(minRGB, input, N*sizeof(double));
    else {
        for (int n=0; n < N; n++)
            if (!udcp)
                minRGB[n] = (input[n] < input[n+N])?((input[n] < input[n+2*N])?(input[n]):(input[n+2*N])):((input[n+N] < input[n+2*N])?(input[n+N]):(input[n+2*N]));
            else
                minRGB[n] = (input[n+N] < input[n+2*N])?(input[n+N]):(input[n+2*N]);
    }

    double *aux = boundaryHandlingImage(minRGB, width, height, sr);
    int widthb = width+2*sr;
#pragma omp parallel for schedule(dynamic)
    for (int j=0; j < height; j++) {
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i < width; i++) {
            double minvalue = aux[(i+sr)+(j+sr)*widthb];
            for (int kj=-sr; (kj <= sr) && (minvalue > 0); kj++)
                for (int ki=-sr; (ki <= sr) && (minvalue > 0); ki++)
                    if (aux[(i+sr+ki)+(j+sr+kj)*widthb] < minvalue) minvalue = aux[(i+sr+ki)+(j+sr+kj)*widthb];

            Jdark[i+j*width] = minvalue;
        }
    }

    delete[] minRGB;
    delete[] aux;

    return Jdark;
}


struct valuesData {
  double value;
  int pixel;
};

//sort data values from highest to lowest
bool sortIncreasing(const valuesData &a, const valuesData &b)
{
    return a.value < b.value;
}

//This function implements Algorithm 1 in accompanying IPOL article
//Input parameters:
//  - dark channel image
//  - input image
//  - width, height and number of channels of input image
//  - p: percentage of pixels used to estimate the ambient light
//  - excludesat: boolen, if TRUE exclude saturated pixels in ambient light estimation
//Output:
//  - ambient light (Ar, Ag, Ab)  (RGB value)
//  - image where pixels used for the computation of the ambient light are marked in red,
//    the rest of pixels keep the values of the input image
double *estimateAmbientLight(double *Jdark, double *input, int width, int height, int nchannels,
                             float p, double &Ar, double &Ag, double &Ab, bool excludesat)
{
    int N = width*height; //total number of pixels
    //fill vector and then sort values (increasing order)
    vector<valuesData> values;
    for (int n=0; n < N; n++) {
        struct valuesData data;
        data.value = Jdark[n];
        data.pixel = n;
        values.push_back(data);
    }
    sort(values.begin(), values.end(), sortIncreasing);

    //mark with Red (black for 1 channel images) pixels used for the computation of the ambient light
    double *Alight = new double[N*nchannels];
    for (int n=0; n < N*nchannels; n++) Alight[n] = input[n];

    //compute average color for brightest pixels in the Jdark image
    Ar = 0;
    Ag = 0;
    Ab = 0;
    int j = 0; //counter for pixels used in the estimation
    int M = floor(p * (float) N /100.0 ); //number of pixels used in the estimation
    for (int i = N-1; (i >= 0) && (j < M); i--) {
        int n = values[i].pixel;
        bool validpixel = true;
        if (excludesat) {
            validpixel = validpixel && (input[n] < 255);
            if (nchannels > 1) {
                validpixel = validpixel && (input[n+N] < 255);
                validpixel = validpixel && (input[n+2*N] < 255);
            }
        }
        if (validpixel) {
            Ar += input[n];
            if (nchannels > 1) {
                Ag += input[n+N];
                Ab += input[n+2*N];
            }
            j++;
            //mark with Red (black for 1 channel images) pixels used for the computation of the ambient light
            Alight[n] = (nchannels > 1)?(255):(0);
            if (nchannels > 1) {
                Alight[n+N] = 0;
                Alight[n+2*N] = 0;
            }
        }
    }

    if (j == 0) {
        //if all pixels are saturated use brightest one
        int n = values[N-1].pixel;
        Ar = input[n];
        Ag = input[n+N];
        Ab = input[n+2*N];
    } else {
        Ar /= (double) j;
        if (nchannels > 1) {
            Ag /= (double) j;
            Ab /= (double) j;
        }
    }

    return Alight;
}


//This function implements line 3 of Algorithm 5 in accompanying IPOL article
//Input parameters:
//  - input image
//  - width, height and number of channels of input image
//  - sr: radius of the square patches used to compute the image of minimum values
//        needed to estimate the transmission map
//  - Ar, Ag, Ab: ambient light (RGB value)
//  - omega: correction factor for transmission map estimation
//  - udcp: boolean flag, if TRUE compute the image of minimum values as minimum over G and B values,
//                        else compute it as minimum over R, G and B
//Output:  transmission map image
double *getTransmissionMap(double *input, int width, int height, int nchannels, int sr,
                           double Ar, double Ag, double Ab, double omega, bool udcp)
{
    int N = width*height;
    double *tmap = new double[N];

    //compute image of minimum values over R, G and B
    double *minRGB = new double[N];
    if (nchannels == 1) {
        for (int n=0; n < N; n++)
            minRGB[n] = input[n]/Ar;
    } else {
        for (int n=0; n < N; n++)
            if (!udcp)
                minRGB[n] = (input[n]/Ar < input[n+N]/Ag)?((input[n]/Ar < input[n+2*N]/Ab)?(input[n]/Ar):(input[n+2*N]/Ab)):((input[n+N]/Ag < input[n+2*N]/Ab)?(input[n+N]/Ag):(input[n+2*N]/Ab));
            else
                minRGB[n] = (input[n+N]/Ag < input[n+2*N]/Ab)?(input[n+N]/Ag):(input[n+2*N]/Ab);
    }


    double *aux = boundaryHandlingImage(minRGB, width, height, sr);
    int widthb = width+2*sr;
#pragma omp parallel for schedule(dynamic)
    for (int j=0; j < height; j++) {
#pragma omp parallel for schedule(dynamic)
        for (int i=0; i < width; i++) {
            double minvalue = aux[(i+sr)+(j+sr)*widthb];
            for (int kj=-sr; (kj <= sr) && (minvalue > 0); kj++)
                for (int ki=-sr; (ki <= sr) && (minvalue > 0); ki++)
                    if (aux[(i+sr+ki)+(j+sr+kj)*widthb] < minvalue) minvalue = aux[(i+sr+ki)+(j+sr+kj)*widthb];

            tmap[i+j*width] = 1- omega * minvalue;
        }
    }

    delete[] minRGB;
    delete[] aux;

    return tmap;
}


//This function implements line 4 of Algorithm 5 in accompanying IPOL article
//Input parameters:
//  - input image
//  - width, height and number of channels of input image
//Output:  intensity image (average of R, G and B values)
double *getIntensity(double *input, int width, int height, int nchannels)
{
    int N = width*height;
    double *I = new double[N];

    if (nchannels == 1) memcpy(I, input, N*sizeof(double));
    else {
        for (int n=0; n < N; n++)
            I[n] = (input[n] + input[n+N] + input[n+2*N]) / 3.0;
    }

    return I;
}


//This function implements Algorithm 2 in accompanying IPOL article
//The code is adapted from the one in IPOL article
//Gabriele Facciolo, Nicolas Limare, and Enric Meinhardt-Llopis,
//Integral Images for Block Matching, Image Processing On Line, 4 (2014), pp. 344–369.
//https://doi.org/10.5201/ipol.2014.57
//The function is based on the function 'computeIntegralImage' from integralImage4c.c
//whose original authors are Gabriele Facciolo and Nicolas Limare
//Input parameters:
//  - input image
//  - width, height input image
//Output:  integral image
double *computeIntegralImage(double *u, int width, int height)
{
   double *ii = new double[width*height];

   // compute the integral image by the recurrence:
   //    s(x,y)  = s(x-1,y) + u(x,y)
   //   ii(x,y) = ii(x,y-1) + s(x,y)
   // with s(0,y) = u(0,y) ; and ii(0,0) = u(0,0)

   // first cell
   ii[0] = u[0];
   // first row
   for (int x = 1; x < width; x++)
      ii[x] = ii[x - 1] + u[x];
   // next rows
   for (int y = 1; y < height; y++) {
      // first cell of the row
      double s = u[width * y];	// current row cumulative sum
      ii[width * y] = ii[width * (y - 1)] + s;
      // next cells
      for (int x = 1; x < width; x++) {
	    s = s + u[x + width * y];
	    ii[x + width * y] = ii[x + width * (y - 1)] + s;
      }
   }

   return ii;
}


//This function implements Algorithm 3 in accompanying IPOL article
//Input parameters:
//  - input image
//  - width, height input image
//  - rr: radius of patch
//Output: image that contains, for each pixel, the average of input values
//        over a patch of size 2*rr+1 x 2*rr+1 centered at the pixel
double *patchAverages(double *I, int width, int height, int rr)
{
    double *Iavg = new double[width*height];

    double *Ib =  boundaryHandlingImage(I, width, height, rr);
    double *II = computeIntegralImage(Ib, width+2*rr, height+2*rr);
    int widthII = width+2*rr;

    int patchsize = (2*rr+1)*(2*rr+1);

#pragma omp parallel for schedule(dynamic)
    for (int j=0; j < height; j++)
        for (int i=0; i < width; i++) {

            int ii = i + rr;
            int jj = j + rr;

            Iavg[i+j*width] = II[(ii+rr)+(jj+rr)*widthII];

            if (ii-rr-1 >= 0)
                Iavg[i+j*width] -= II[(ii-rr-1)+(jj+rr)*widthII];

            if (jj-rr-1 >= 0)
                Iavg[i+j*width] -= II[(ii+rr)+(jj-rr-1)*widthII];

            if ((ii-rr-1 >= 0) && (jj-rr-1 >= 0))
                Iavg[i+j*width] += II[(ii-rr-1)+(jj-rr-1)*widthII];

            Iavg[i+j*width] /= (double) patchsize;
        }

    delete[] Ib;
    delete[] II;

    return Iavg;
}


//This function implements Algorithm 4 (guided filter) in accompanying IPOL article
//Input parameters:
//  - original transmission map
//  - input image (guiding image)
//  - width, height input image
//  - rr: radius of the square patches used in the guided filter
//  - eps: regularization parameter of the guided filter
//Output: refined transmission map
double *refineTransmissionMap(double *tmap, double *I, int width, int height, int rr, double eps)
{
    int N = width*height;

    //normalize I to range [0, 1]
    for (int n=0; n < N; n++) I[n] /= 255.0; //warning: the values in array I are modified!

    double *mu_I = patchAverages(I, width, height, rr);

    double *I2 = new double[N];
    for (int n=0; n < N; n++) I2[n] = I[n] * I[n];
    double *mu_I2 = patchAverages(I2, width, height, rr);

    double *var_I = new double[N];
    for (int n=0; n < N; n++) var_I[n] = mu_I2[n] - mu_I[n] * mu_I[n];

    double *mu_t = patchAverages(tmap, width, height, rr);

    double *It = new double[N];
    for (int n=0; n < N; n++) It[n] = I[n] * tmap[n];
    double *mu_It = patchAverages(It, width, height, rr);

    double *covar_It = new double[N];
    for (int n=0; n < N; n++) covar_It[n] = mu_It[n] - mu_I[n] * mu_t[n];

    double *a = new double[N];
    for (int n=0; n < N; n++) a[n] = covar_It[n] / (var_I[n] + eps);

    double *b = new double[N];
    for (int n=0; n < N; n++) b[n] = mu_t[n] - a[n] * mu_I[n];

    double *mu_a = patchAverages(a, width, height, rr);
    double *mu_b = patchAverages(b, width, height, rr);

    double *tmapf = new double[N];
    for (int n=0; n < N; n++) tmapf[n] = mu_a[n] * I[n] + mu_b[n];


    delete[] mu_I;
    delete[] I2;
    delete[] mu_I2;
    delete[] var_I;
    delete[] mu_t;
    delete[] It;
    delete[] mu_It;
    delete[] covar_It;
    delete[] a;
    delete[] b;
    delete[] mu_a;
    delete[] mu_b;

    return tmapf;
}



//This function implements line 6 of Algorithm 5 in accompanying IPOL article
//Input parameters:
//  - input image
//  - width, height and number of channels of input image
//  - Ar, Ag, Ab: ambient light (RGB value)
//  - refined transmission map
//  - t0: minimum allowed value of the transmission map for recovery of dehazed result
//Output: dehazed image
double *getRadiance(double *input, int width, int height, int nchannels,
                    double Ar, double Ag, double Ab, double *tmapf, double t0)
{
    int N = width*height;
    double *dehazed = new double[N*nchannels];

    double A, t;
    for (int n=0; n < N; n++) {
        t = (tmapf[n] > t0)?(tmapf[n]):(t0);
        for (int k=0; k < nchannels; k++) {
            A = (k == 0)?(Ar):((k == 1)?(Ag):(Ab));
            dehazed[n+k*N] = (input[n+k*N]/255.0 - A) / t + A;
        }
    }

    return dehazed;
}

//Main algorithm: this function implements Algorithm 5 in accompanying IPOL article
int main(int argc,  char **argv)
{
    //Get input parameters
    std::vector <OptStruct *> options;
    OptStruct os = {(char *) "s:", 0, (char *) "7", NULL, (char *) "Radius of the square patches used to compute the dark channel"};
    options.push_back(&os);
    OptStruct op = {(char *) "p:", 0, (char *) "0.1", NULL, (char *) "Percentage of pixels used to estimate the ambient light"};
    options.push_back(&op);
    OptStruct ow = {(char *) "w:", 0, (char *) "0.95", NULL, (char *) "Correction factor for transmission map estimation"};
    options.push_back(&ow);
    OptStruct og = {(char *) "r:", 0, (char *) "30", NULL, (char *) "Radius of the square patches used in the guided filter"};
    options.push_back(&og);
    OptStruct oe = {(char *) "e:", 0, (char *) "0.0001", NULL, (char *) "Regularization parameter of the guided filter"};
    options.push_back(&oe);
    OptStruct ot = {(char *) "t:", 0, (char *) "0.1", NULL, (char *) "Minimum allowed value of the transmission map for recovery of dehazed result"};
    options.push_back(&ot);

    //Optional parameters
    OptStruct oa = 	{(char *) "a:", 0, NULL, NULL, (char *) "ambient light image (optional output)"};
    options.push_back(&oa);
    OptStruct od = 	{(char *) "d:", 0, NULL, NULL, (char *) "dark channel image (optional output)"};
    options.push_back(&od);
    OptStruct om = 	{(char *) "m:", 0, NULL, NULL, (char *) "transmission map image (optional output)"};
    options.push_back(&om);
    OptStruct of = 	{(char *) "f:", 0, NULL, NULL, (char *) "refined transmission map image (optional output)"};
    options.push_back(&of);
    OptStruct ou = {(char *) "u", 0,  NULL, NULL, (char *) "Underwater DCP"};
    options.push_back(&ou);
    OptStruct ox = {(char *) "x", 0,  NULL, NULL, (char *) "Exclude saturated pixels in ambient light estimation"};
    options.push_back(&ox);

    //Input and output images
    std::vector<ParStruct *> pparameters;
    ParStruct pinput = {(char *) "input", NULL, (char *) "input image"};
    pparameters.push_back(&pinput);
    ParStruct poutput = {(char *) "output", NULL, (char *) "output image"};
    pparameters.push_back(&poutput);

    if (!parsecmdline((char *) "dehazeDCP", (char *) "Dehazing based on Dark Channel Prior",
                      argc, argv, options, pparameters))
        return EXIT_FAILURE;


    int sr = atoi(os.value); //Radius of the square patches used to compute the dark channel
    int rr = atoi(og.value); //Radius of the square patches used in the guided filter
    float p = atof(op.value); //Percentage of pixels used to estimate the ambient light
    double omega = (double) atof(ow.value); //Correction factor for transmission map estimation
    double eps = (double) atof(oe.value); //Regularization parameter of the guided filter
    double t0 = (double) atof(ot.value); //Minimum allowed value of the transmission map for recovery of dehazed result

    const char *namealight=NULL;
    if (oa.flag) namealight = oa.value;
    const char *namedark=NULL;
    if (od.flag) namedark = od.value;
    const char *nametmap=NULL;
    if (om.flag) nametmap = om.value;
    const char *nametmapf=NULL;
    if (of.flag) nametmapf = of.value;

    bool udcp = ou.flag;
    bool excludesat = ox.flag;

    char *namein = pinput.value;
    char *nameout = poutput.value;

    //Read input data
    float *dataIn;
    int width, height, nchannels;
    dataIn=iio_read_image_float_split(namein, &width, &height, &nchannels);

    printf("%i x %i   %i channels\n", width, height, nchannels);

    //do not use alpha channel
    if(nchannels >= 3) nchannels=3;
    else nchannels=1;

    //convert input data to double
    double *input = new double[width*height*nchannels];
    for (int n=0; n < width*height*nchannels; n++) input[n] = (double) dataIn[n];


    //Estimate Dark Channel
    double *Jdark = getDarkChannel(input, width, height, nchannels, sr, udcp);

    //Estimate ambient light
    double Ar, Ag, Ab;
    double *Alight = estimateAmbientLight(Jdark, input, width, height, nchannels, p, Ar, Ag, Ab, excludesat);

    printf("Estimated ambient light: (%2.2f, %2.2f, %2.2f)\n", Ar, Ag, Ab);
    if (excludesat) printf("Excluded saturated pixels\n");

    //Compute transmission map
    double *tmap = getTransmissionMap(input, width, height, nchannels, sr, Ar, Ag, Ab, omega, udcp);

    //Compute Intensity
    double *I = getIntensity(input, width, height, nchannels);

    //Refine transmission map (use guided filter)
    double *tmapf = refineTransmissionMap(tmap, I, width, height, rr, eps);

    //Recover scece radiance
    double *dehazed = getRadiance(input, width, height, nchannels, Ar/255.0, Ag/255.0, Ab/255.0, tmapf, t0);


    //output files
    //reuse dataIn data structure
    for (int n=0; n < width*height*nchannels; n++) dataIn[n] = (float) clip(dehazed[n]*255.0, 0, 255);
    iio_save_image_float_split(nameout, dataIn, width, height, nchannels);

    if (namealight) {
        for (int n=0; n < width*height*nchannels; n++) dataIn[n] = (float) Alight[n];
        iio_save_image_float_split((char *) namealight, dataIn, width, height, nchannels);
    }

    if (namedark) {
        for (int n=0; n < width*height; n++) dataIn[n] = (float) Jdark[n];
        iio_save_image_float_split((char *) namedark, dataIn, width, height, 1);
    }

    if (nametmap) {
        for (int n=0; n < width*height; n++) dataIn[n] = (float) tmap[n]*255.0;
        iio_save_image_float_split((char *) nametmap, dataIn, width, height, 1);
    }

    if (nametmapf) {
        for (int n=0; n < width*height; n++) dataIn[n] = (float) tmapf[n]*255.0;
        iio_save_image_float_split((char *) nametmapf, dataIn, width, height, 1);
    }


    //clean up
    delete[] input;
    delete[] Jdark;
    delete[] tmap;
    delete[] tmapf;
    delete[] Alight;
    delete[] dehazed;
    delete[] I;
    free(dataIn);

    return EXIT_SUCCESS;

}




