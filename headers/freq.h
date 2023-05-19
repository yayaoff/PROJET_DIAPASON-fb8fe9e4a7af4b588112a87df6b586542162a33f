#include "eigen.h"
#include <math.h>
#include <stdlib.h>
#include "lu.h"
#include "elasticity.h"
// #include "../gmsh-4.11.1-Windows64-sdk/include/gmshc.h"
#include "design.h"
#include "RCM.h"

double *compute_freq(int tag_vis,param_t *param);

double *compute_band_freq(int tag_vis,param_t *param);

double *compute_band_freq_half(int tag_vis,param_t *param);