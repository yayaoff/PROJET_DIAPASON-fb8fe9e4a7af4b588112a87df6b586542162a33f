#include "../headers/opti.h"

int main(int argc, char *argv[]){

    if (argc < 6){
    printf("Usage: \n"
			"./opti <k> <out> <r1> <r2> <e> <l> <out>\n" 
			"---------------------------- \n\n"
      "- k\n "
      "- out\n "
			"- r1\n "
			"- r2\n "
      "- e\n "
      "- l\n "
      "\n");
		return -1;
  } 
  param_t *params = malloc(sizeof(param_t));
  params->k = atoi(argv[1]); params->argc = argc; params->argv = argv; 
  params->r1 = atof(argv[3]) ; params->r2 = atof(argv[4]); params->e= atof(argv[5]);params->l = atof(argv[6]);
  params->meshSizeFactor=0.3;params->filename=argv[2];
  
  double *freqs = compute_band_freq(0,params);

  return 0;
}