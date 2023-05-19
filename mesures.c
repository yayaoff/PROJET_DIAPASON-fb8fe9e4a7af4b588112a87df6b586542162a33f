#include <stdlib.h>
#include <string.h>
// #include "../gmsh-sdk/include/gmshc.h"
#include "../gmsh-4.11.1-Windows64-sdk/include/gmshc.h"
#include "headers/matrix.h"
#include "headers/elasticity.h"
#include "math.h"
#include "headers/lu.h"
#include "headers/design.h"
#include "headers/opti.h"
// #include "headers/freq.h"
#include <time.h>

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

void write_csv(double time, double mesh_size, char* filename){
    FILE* file = fopen(filename, "a");
    if (file == NULL) {
        printf("Error opening the file.\n");
        return;
    }

    fprintf(file, "%f,%f\n", mesh_size, time);

    fclose(file);
}

int main (int argc, char *argv[]) {

  if (argc < 3){
    printf("Usage: \n"
			"./project <k> <out> <half>\n" 
			"---------------------------- \n\n"
			"- k is the number of frequencies to compute. \n "
			"- out is the output file to write the frequencies. \n "
      "- <half> is 1 if half designed, 0 otherwise"
      "\n");
		return -1;
  } 

  param_t *param_init = malloc(sizeof(param_t));
  param_init->argc = argc; param_init->argv = argv; param_init->k = atoi(argv[1]);
  param_init->r1=0.006415610882232378; param_init->r2=0.015; param_init->e=0.058552936021614585; param_init->l=0.05740262518415047; param_init->meshSizeFactor=0.23;param_init->filename=NULL;
  // clock_t start, end;
  // double cpu_time_used;

  //   start = clock();  // Record the start time
       compute_freq(0,param_init);
  //   end = clock();  // Record the end time
  //   cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //   write_csv(cpu_time_used,param_init->meshSizeFactor,"plots/matrix.csv");

    // start = clock();  // Record the start time
    // compute_freq_half(0,param_init);
    // end = clock();  // Record the end time
    // cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    // write_csv(cpu_time_used,param_init->meshSizeFactor,"plots/matrix_half.csv");

    // start = clock();  // Record the start time
    // compute_band_freq(0,param_init);
    // end = clock();  // Record the end time
    // cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    // write_csv(cpu_time_used,param_init->meshSizeFactor,"plots/matrix_band.csv");

    // start = clock();  // Record the start time
    // compute_band_freq_half(0,param_init);
    // end = clock();  // Record the end time
    // cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    // write_csv(cpu_time_used,param_init->meshSizeFactor,"plots/matrix_band_half.csv");

    return 0;
}