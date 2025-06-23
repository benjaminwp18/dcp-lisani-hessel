### Source code for IPOL article ["Dehazing with Dark Channel Prior: Analysis and Implementation"](https://ipol.im)


Version 1.1 - April 20, 2024
by Jose-Luis Lisani <joseluis.lisani@uib.es>


### ASSOCIATED PUBLICATION

This code is associated to the following IPOL publication:

J.L. Lisani, Charles Hessel,
Dehazing with Dark Channel Prior: Analysis and Implementation
Image Processing On Line (IPOL), https://www.ipol.im/pub/art/2024/530/

An online demo that uses the code is available at IPOL: https://ipolcore.ipol.im/demo/clientApp/demo.html?id=530


### ORGANIZATION OF THE CODE

The code (in srcIPOL/src folder) is organized as follows:

- a folder 'library' containing the following: 
	 - libImageFormats.h, libImageFormats.c: auxiliary functions for 
   images input/output
     - parser.h, parser.cpp: functions to manage the input parameters of the
   main function
- dehazeDCP.cpp: the main function
- makefile: for compilation
- agpl-3.0.txt: GNU Affero General Public License 

Moreover, two test images are provided (in srcIPOL/example_images):
- airfield.png (credit: NASA http://dragon.larc.nasa.gov/retinex/)
- underwater.png (credit: IMEDEA)

The results of the examples of use (at the end of this document) are stored, 
as reference, in srcIPOL/results_examples.

The algorithm described in the associated paper is implemented 
in dehazeDCP.cpp. This code has been peer-reviewed by IPOL.

### IMPORTANT NOTES

- The code has been compiled and tested using a Ubuntu Linux environment
- On Ubuntu you may need to install the following libraries: libjpeg, libpng, libtiff
````{verbatim}
sudo apt-get install build-essential libjpeg8-dev libpng-dev libtiff-dev
````
 

### COMPILATION

1) Decompress code and access to folder:
````bash
tar xvzf srcIPOL.tgz
cd srcIPOL/src
````

3) Compilation:
````bash
make OMP=1
````

The executable file obtained after compilation is 'dehazeDCP'.



### USAGE

````bash
usage: dehazeDCP [-s s] [-p p] [-w w] [-r r] [-e e] [-t t] [-a a] [-d d] [-m m] [-f f] [-u] [-x]  input  output 
	-s  s	 Radius of the square patches used to compute the dark channel (Default: 7)
	-p  p	 Percentage of pixels used to estimate the ambient light (Default: 0.1)
	-w  w	 Correction factor for transmission map estimation (Default: 0.95)
	-r  r	 Radius of the square patches used in the guided filter (Default: 30)
	-e  e	 Regularization parameter of the guided filter (Default: 0.0001)
	-t  t	 Minimum allowed value of the transmission map for recovery of dehazed result (Default: 0.1)
	-a  a	 ambient light image (optional output) 
	-d  d	 dark channel image (optional output) 
	-m  m	 transmission map image (optional output) 
	-f  f	 refined transmission map image (optional output) 
	-u	 Underwater DCP 
	-x	 Exclude saturated pixels in ambient light estimation 
	input	 input image
	output	 output image
````




### EXAMPLES OF USE

The following commands apply dehazeDCP to the provided image (run from srcIPOL folder):

1) Dark Channel Prior dehazing (default parameters, only output dehazed image)

````bash
./src/dehazeDCP ./example_images/airfield.png airfied_output.png
````

2) Dark Channel Prior dehazing (default parameters, output dehazed image and auxiliary results)

````bash
./src/dehazeDCP -a airfied_alight.png -d airfied_dark.png -m airfied_tmap.png -f airfied_tmapf.png ./example_images/airfield.png airfied_output.png 
````

3) Dark Channel Prior dehazing (parameters different from default, only output dehazed image)

````bash
./src/dehazeDCP -s 10 -p 0.2 -w 0.9 -r 50 -e 0.01 -t 0.2 ./example_images/airfield.png airfied_output2.png
````

4) Underwater Dark Channel Prior dehazing (default parameters, output dehazed image and auxiliary results)

````bash
./src/dehazeDCP -a underwater_alight.png -d underwater_dark.png -m underwater_tmap.png -f underwater_tmapf.png -u ./example_images/underwater.png underwater_output.png 
````



