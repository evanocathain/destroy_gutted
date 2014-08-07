#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"

void print_usage(void);

int get_args(int argc, char** argv, char filename[], int* ascii, int* sigproc, int* nsamps, float* spthresh, int* nsmax, int* subzero, double* dm, int* boxcar_option)
{
  if (argc==1){
    print_usage();
  }

  fprintf(stdout,"DESTROY - version 0.1\n\n"); // Program launch ...

  /* Initialise the default values */
  *ascii=0; 
  *sigproc=1;
  *spthresh=4.00;
  *dm=0.0;
  *nsmax=10;
  *boxcar_option=1;

  /* Get values from header */
  /* ADD THIS LATER */

  /* Cycle through all the command line arguments */
  int i=1;
  while (i<argc){
    if (strings_equal(argv[i],"-o")){
      i++;
    }else if (strings_equal(argv[i],"-ascii")){
      *ascii=1;
      *sigproc=0;
    }else if (strings_equal(argv[i],"-box")){
      i++;
      *boxcar_option=atoi(argv[i]);
    }else if (strings_equal(argv[i],"-spthresh")){
      i++;
      *spthresh=atof(argv[i]);
    }else if (strings_equal(argv[i],"-dm")){
      i++;
      *dm=atof(argv[i]);
    }else if (strings_equal(argv[i],"-n")){
      i++;
      *nsamps=atoi(argv[i]);
    }else if (strings_equal(argv[i],"-nsmax")){
      i++;
      *nsmax=atoi(argv[i]);
    }else if (strings_equal(argv[i],"-h")){
      print_usage();
    }else if (strings_equal(argv[i],"--help")){
      print_usage();
    }
    i++;
  }
  strcpy(filename,argv[argc-1]);

  return(0);
}

void print_usage(void)
{
  fprintf(stdout,"\ndestroy - 'seeks' out transient signals in a noisey time series\n\n");
  fprintf(stdout,"usage: destroy -{options} {filename}\n\n");
  fprintf(stdout,"options:\n\n");
  fprintf(stdout,"-ascii          - set input file format to ascii (def=sigproc)\n");
  fprintf(stdout,"-dm val         - set dm to val (def=0.0 or from sigproc header)\n");
  fprintf(stdout,"Single Pulse Search Options\n");
  fprintf(stdout,"-box val        - set boxcar options: 0=decimate, pow of 2 steps; 1=convolve, pow of 2 steps; 2=convolve, steps of 1 (def=1)\n");
  fprintf(stdout,"-spthresh val   - set search threshold to val sigma (def=4.0)\n");
  fprintf(stdout,"-n val          - set the number of samples to read in (def=all)\n");
  fprintf(stdout,"-nsmax val      - set the highest boxcar size to search to be 2**nsmax (def=10)\n");
  fprintf(stdout,"-maura          - output file with highest S/N detections only (def=don't)[DOESN'T WORK YET!]\n");
  fprintf(stdout,"\n-h,--help       - prints this usage message\n");
  fprintf(stdout,"\n");

  exit(0);
}
