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

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

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

  int tag_half = atoi(argv[3]);
  printf("tag half = : %d\n",tag_half);

  // We want fabs(freq-f_target) <= EPS

  // Parameters to optimize : r1, r2, e, l
  // Ici, on optimise le paramètre r1 et on fixe les autres arbitrairement-> Tester pour les 4 et voir
  // pour quel paramètre on converge le plus vite ?
  // Rmq : attention aux contraintes : r1 < r2
  param_t *param_init = malloc(sizeof(param_t));
  param_init->argc = argc; param_init->argv = argv; param_init->k = atoi(argv[1]);
  param_init->r1=0.006415610882232378; param_init->r2=0.015; param_init->e=0.058552936021614585; param_init->l=0.05740262518415047; param_init->meshSizeFactor=0.3;param_init->filename=NULL;
  // Compute frequencies 
  // If half designed:
  if (tag_half==1){
    double* new_frequencies = compute_band_freq_half(2,param_init);
    printf("Freq final = %f\n",new_frequencies[0]);
    
  }
  else{
      double* new_frequencies = compute_band_freq(2,param_init);
      printf("Freq final = %f\n",new_frequencies[0]);
  }
  return 0;
}