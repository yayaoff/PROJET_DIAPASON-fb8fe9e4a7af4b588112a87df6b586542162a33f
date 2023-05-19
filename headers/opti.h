#include "freq.h"

#define TOL 5e-6
#define f_target 1244.51    //frequence for D5# [Hz] 

double *compute_freq(int tag_vis,param_t *param);
double *compute_band_freq(int tag_vis,param_t *param);
double f(param_t* param, double param_test);
double bissection_method(double low_bound, double up_bound,param_t *param);

double bissection_half_method(double low_bound, double up_bound,param_t *param);
double f_half(param_t* param, double param_test);