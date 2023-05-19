#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "../gmsh-sdk/include/gmshc.h"
#include "gmsh-4.11.1-Windows64-sdk/include/gmshc.h"
#include "headers/matrix.h"
#include "headers/elasticity.h"
#include "math.h"
#include "headers/lu.h"
#include "headers/design.h"
#include "headers/eigen.h"
#include "headers/opti.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

// int k;
// double r1; double r2; double e; double l; double meshSizeFactor; char * filename;

int main (int argc, char *argv[]) {

  if (argc < 2){
    printf("Usage: \n"
			"./project <k> <out>\n" 
			"---------------------------- \n\n"
			"- k is the number of frequencies to compute. \n "
			"- out is the output file to write the frequencies. \n "
      "\n");
		return -1;
  } 
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);
  designTuningFork(0.006415610882232378,0.015,0.058552936021614585,0.05740262518415047,0.3,NULL);
  return 0;
}