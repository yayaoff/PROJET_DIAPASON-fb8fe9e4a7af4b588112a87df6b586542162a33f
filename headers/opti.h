#include "eigen.h"
#include <math.h>
#include <stdlib.h>
#include "lu.h"
#include "elasticity.h"
#include "../../gmsh-sdk/include/gmshc.h"
#include "design.h"
#include "RCM.h"

#define TOL 1e-10
#define f_target 1244.51    //frequence for D5# [Hz] 

typedef struct{
    int argc;
    int k;
    double r1;
    double r2;
    double e;
    double l; 
    double meshSizeFactor;
    char *filename;
    char **argv;
}param_t;

double *compute_freq(int tag_vis,param_t *param);
double *compute_band_freq(int tag_vis,param_t *param);
double f(param_t* param, double param_test);
double bissection_method(double low_bound, double up_bound,param_t *param);