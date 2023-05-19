#include "../headers/opti.h"
// Méthode basique : méthode de bissection 
// Défaut : Optimisation univariée, or plusieurs paramètres !!
// On fixe 3 paramètres sur 4 et on effectue la méthode sur le dernier paramètre
// Amélioration potentielle : Utiliser une méthode d'optimisation multivariée (descente de gradient?) 

double bissection_method(double low_bound, double up_bound,param_t *param){
    printf("\n------------------------------------------------------\n");
    printf(" Bissection method !\n");
    printf("------------------------------------------------------\n");
    

    double freq_mid, mid;
    double act = param->l;
    param->l = low_bound;
    double freq_low = compute_band_freq(0,param)[0];
    param->l = up_bound;
    double freq_up = compute_band_freq(0,param)[0];

    printf("freq_low computed : %f\n",freq_low);
     printf("freq_up computed : %f\n",freq_up);

    if((freq_low - f_target) * (freq_up - f_target)>0){
      printf("Non opposite signs\n");
      return -1;
    }

    while (fabs(up_bound - low_bound) >  TOL) {
      mid = fabs(up_bound + low_bound) / 2.0;
      param-> l = mid;
      freq_mid = compute_band_freq(0,param)[0];
      if((freq_low - f_target)*(freq_mid - f_target) < 0){
        up_bound = mid;
      }
      else{
        low_bound = mid;
        param->l = low_bound;
        freq_low = compute_band_freq(0,param)[0];
      }
      param->l = act;
    }
    return mid;
}

double bissection_half_method(double low_bound, double up_bound,param_t *param){
    printf("\n------------------------------------------------------\n");
    printf(" Bissection method !\n");
    printf("------------------------------------------------------\n");

    double freq_mid, mid;
    double act = param->l;
    param->l = low_bound;
    double freq_low = compute_band_freq_half(0,param)[0];
    param->l = up_bound;
    double freq_up = compute_band_freq_half(0,param)[0];

    if((freq_low - f_target) * (freq_up - f_target)>0){
      printf("Non opposite signs\n");
      return -1;
    }

    while (fabs(up_bound - low_bound) >  TOL) {
      mid = fabs(up_bound + low_bound) / 2.0;
      param-> l = mid;
      freq_mid = compute_band_freq_half(0,param)[0];
      if((freq_low - f_target)*(freq_mid - f_target) < 0){
        up_bound = mid;
      }
      else{
        low_bound = mid;
        param->l = low_bound;
        freq_low = compute_band_freq_half(0,param)[0];
      }
      param->l = act;
    }
    return mid;
}