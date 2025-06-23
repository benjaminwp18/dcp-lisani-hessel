/*
 *  parser.cpp
 *
 *
 *  Created by Antoni Bueades Capó.
 * 
 *  Copyright 2012 IPOL Image Processing On Line http://www.ipol.im/
 *
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

 


#include "parser.h"



void  printusage(char *pname, char *gp,std::vector<OptStruct*>  &opt,
				 std::vector<ParStruct*>  &par)
{
        
        //int nopt = opt.size();
        int npar = par.size();
        
        ///// USAGE
        printf("\nusage: %s ", pname);
        for(int i=0; i < (int) strlen(gp); i++)
                if (gp[i] != ':')
                {
                        printf("[-%c",gp[i]);
                        
                        if (i+1 < (int) strlen(gp) && gp[i+1] ==  ':') printf(" %c] ", gp[i]);
                        else printf("] ");
                        
                }
        
        for(int i=0; i < npar; i++)
                printf(" %s ", par[i]->name);
        
        printf("\n");
        //// PARAMETERS
        
        int j=0;
        for(int i=0; i < (int) strlen(gp); i++)
                if (gp[i] != ':')
                {
                        printf("\t-%c",gp[i]);
                        
                        if (i+1 < (int) strlen(gp) && gp[i+1] ==  ':') {
                                
                                printf("  %c\t %s ", gp[i], opt[j]->comment);
                                if (opt[j]->defvalue != NULL) printf("(Default: %s)",opt[j]->defvalue);
                                
                                printf("\n");
                                
                        }               
                        else printf("\t %s \n", opt[j]->comment);
                        
                        j++;
                }
        
        
        for(int i=0; i < npar; i++)
        {
                printf("\t%s",par[i]->name);
                printf("\t %s\n", par[i]->comment);
        }
        
        
}




int parsecmdline(char *pname, char *function, int argc, char **argv,
                 std::vector <OptStruct*> & opt, std::vector <ParStruct*> & par)
{
        int nopt = opt.size();
        int npar = par.size();
        
        
        char *gp = new char[2*nopt+1];
        gp[0]='\0';
        
        
        for(int i=0; i < nopt; i++) {
			opt[i]->flag = 0; 
			opt[i]->value=NULL; 
			strcat(gp, opt[i]->gp);
		}
        for(int i=0; i < npar; i++) { par[i]->value = NULL;}
        
        opterr = 0;     // No messages by getopt
        
        int c;
        while ((c = getopt (argc, argv, gp)) != -1)
        {
                
                int j=0;
                for(unsigned int i=0; i < strlen(gp); i++)
                        if (c == gp[i])
                        {
                                
                                opt[j]->flag = 1;
                                opt[j]->value = optarg;
                                break;
                                
                        } else if (gp[i] != ':') j++;
                
                
                
                if (c == '?')
                {       
                        
                        unsigned int i = 0;
                        for(i=0; i < strlen(gp); i++)
                                if (optopt == gp[i])
                                {
                                        printf("\n%s: %s\n", pname, function);
                                        fprintf (stderr, "\nerror: option -%c requires an argument.\n", optopt);
                                        break;  
                                }
                        
                        if (i == strlen(gp)) {  printf("\n%s: %s\n", pname, function);
                                fprintf (stderr, "\nerror: unknown option `-%c'.\n", optopt);
                        }
                        
                        printusage(pname, gp,  opt,  par);
                        return 0;
                        
                }
                
        }
        
        
        //// Setting default values for non selected options
        for(int j=0; j < nopt; j++)
                if (opt[j]->flag == 0 && opt[j]->defvalue != NULL)
					opt[j]->value =  opt[j]->defvalue;
        
        
        if (argc - optind != npar) {
                printf("\n%s: %s\n", pname, function);
                fprintf (stderr, "\nerror: incorrect number of parameters\n");
                printusage(pname, gp, opt,par);
                return 0;
        }
        int i=0;
        for (int index = optind; index < argc ; index++, i++){
                par[i]->value = argv[index];
        }
        
        return 1;
        
        
}

