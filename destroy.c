#include <stdio.h>
#include <string.h>
#include "functions.h"

int main(int argc, char* argv[])
{
  int check=0, nsamps=0, ascii, sigproc, subzero, nsmax=10, boxcar_option;
  char filename[100];
  float spthresh, *series;
  double dm=0.0, tstart, tsamp, fch1, foff, refdm, raj, decj;

  check=get_args(argc,argv,filename,&ascii,&sigproc,&nsamps,&spthresh,&nsmax,&subzero,&dm,&boxcar_option); // get cmd line args
  series=read_data(filename,series,&nsamps,sigproc,ascii,&dm,&raj,&decj,&tstart,&tsamp,&fch1,&foff);       // read in the data
  sp_search(series,spthresh,nsmax,nsamps,dm,tsamp,boxcar_option);                                          // search the data

  fprintf(stdout,"Destruction ended ...\n");
  return(0);
}
