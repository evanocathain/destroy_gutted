#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
//#include <limits.h> USED TO CHECK WHAT INT_MAX was during development

void decimate(int* npoints, float* series, int dec_fac);
void find_mean_sig(int nsamps, float* series, float* mean, float* sig, float thresh);
void convolve_boxcar(int* npoints, int nsamps, float* series, float* series_conv, int boxcar_size);

void sp_search(float* series, float thresh, int nsmax, int nsamps, double dm, double tsamp, int boxcar_option)
{

  int i, j=0, ns=0, npoints, npulses=0, phaseok;
  float mean, sig, snr;
  FILE *pulses, *hst;

  fprintf(stdout,"Single Pulse Search:\tThreshold = %.2f sigma\n",thresh);
 
  /* Open output files ... */
  pulses=fopen("pulses.pls","a");
  hst=fopen("pulses.hst","a");
  /* Write header line - ADD LATER */

  /* Search the time series once */ 
  find_mean_sig(nsamps,series,&mean,&sig,thresh); 
  for(i=0; i<nsamps; i++){
    snr=(series[i]-mean)/sig;
    if (snr>=thresh){
      npulses++;
      fprintf(pulses,"%.4f 1 %.1f %.2f\n",dm,i+0.5,snr);                      // Output: DM, boxcar_size, pulse time, S/N
    }
  }
  npoints=nsamps;

  /* Now loop through and search time series for different size box-cars */
  int boxcar=0;
  if (boxcar_option==0){
    // FASTEST WAY, BUT HAS SQRT(2) ISSUE IN S/N CALCS
    for(ns=1; ns<=nsmax; ns++){                                               // Boxcars 2, 4, 8, 16, 32, ...
      boxcar = (int)pow(2,ns);
      decimate(&npoints,series,2);                                            // decimates by factor of 2, saves _in place_
      find_mean_sig(npoints,series,&mean,&sig,thresh);
      //printf("%d %f %f %d\n", npoints, mean, sig, boxcar);
      for(i=0; i<npoints; i++){
	snr=(series[i]-mean)/sig;
	if (snr>=thresh){
	  npulses++;
	  fprintf(pulses,"%.4f %d %d %.2f\n",dm,boxcar,i*boxcar+boxcar/2,snr);// Output: DM, boxcar_size, pulse time, S/N
	}
      }
    }
  }else if (boxcar_option==1){
    // SLOWER WAY, BUT CORRECT FOR POWERS OF 2
    float *series_conv;
    find_mean_sig(npoints,series,&mean,&sig,thresh);
    series_conv=(float*)malloc(nsamps*sizeof(float));                         // CREATE ARRAY 
    for(ns=1; ns<=nsmax; ns++){                                               // Boxcars 2, 4, 8, 16, 32, ...
      boxcar = (int)pow(2,ns);
      convolve_boxcar(&npoints,nsamps,series,series_conv,boxcar);
      for(i=0; i<npoints; i++){
	snr=(series_conv[i]-mean)/(sig*sqrt(boxcar));
	if (snr>=thresh){
	  npulses++;
	  fprintf(pulses,"%.4f %d %d %.2f\n",dm,boxcar,i+boxcar/2,snr);       // Output: DM, boxcar_size, pulse time, S/N
	}
      }
    }
  }else if (boxcar_option==2){  
    // SLOWEST WAY, BUT CORRECT FOR ALL BOXCARS
    float *series_conv;
    find_mean_sig(npoints,series,&mean,&sig,thresh);
    series_conv=(float*)malloc(nsamps*sizeof(float));                         // CREATE ARRAY 
    for(ns=1; ns<=(int)pow(2,nsmax); ns++){                                   // Boxcars 2, 3, 4, 5, 6, 7, 8, 9, ...
      boxcar = ns;
      convolve_boxcar(&npoints,nsamps,series,series_conv,boxcar);
      for(i=0; i<npoints; i++){
	snr=(series_conv[i]-mean)/(sig*sqrt(boxcar));
	if (snr>=thresh){
	  npulses++;
	  fprintf(pulses,"%.4f %d %.1f %.2f\n",dm,boxcar,i+boxcar*0.5,snr);   // Output: DM, boxcar_size, pulse_time, S/N
	}
      }
    }
  }else {
    fprintf(stderr,"Unrecognised value %d for -box option, use 0, 1 or 2\n",boxcar_option);
  }
  fprintf(hst,"%.4f %d\n",dm,npulses);

  /* Close output files ... */
  fclose(pulses);
  fclose(hst);

  /* if -maura option given sort the output pls file to only keep the
     highest S/N event for pulses detected at more than width. Write
     this to pulses.maura.pls or something like that  */
  //  system("sort -n -k 3 pulses.pls USE SORT AND UNIQ"); 

  return;
}

/* Function Definitions */
void find_mean_sig(int nsamps, float* series, float* mean, float* sig, float thresh)
{
  int i,j=0;
  float s,ss;
  s=ss=0.0;
  /* Find mean and rms */
  for(i=0; i<nsamps; i++){
    s+=series[i];
    ss+=(series[i])*(series[i]);
    j++;
  }
  *mean=s/(float)j;
  *sig=sqrt((ss/(float)j)-((*mean)*(*mean)));

  /* Scan through again to remove bright spikes */
  s=ss=0.0;
  j=0;
  for(i=0; i<nsamps; i++){
    if((series[i]-(*mean))/(*sig)<3.0){
      s+=series[i];
      ss+=(series[i])*(series[i]);
      j++;
    }
  }
  *mean=s/(float)j;
  *sig=sqrt((ss/(float)j)-((*mean)*(*mean)));
  return;
}

void decimate(int* npoints, float* series, int dec_fac)
{
  int i, j, n, n_decimated;
  n = *npoints;
  n_decimated = (n/dec_fac);   // integer division, lose n%dec_fac samples

  /* The dec_fac == 2 case, should make it general ... */
  for(i=0; i<n_decimated; i++){
    series[i] = series[2*i] + series[2*i+1];
  }
    
  *npoints=n_decimated;
  return;
}

void convolve_boxcar(int* npoints, int nsamps, float* series, float* series_conv, int boxcar_size)
{
  int i;
  float run_sum=0.0;

  *npoints = nsamps - (boxcar_size - 1);     // you lose (boxcar_size - 1)/2 samples after convolution

  for (i=0; i<boxcar_size; i++){
    run_sum += series[i];
  }
  for (i=boxcar_size; i<*npoints; i++){
    series_conv[i-boxcar_size] = run_sum;
    run_sum += series[i];
    run_sum -= series[i-boxcar_size];
  }
  
}
