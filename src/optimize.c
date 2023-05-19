#include "../headers/opti.h"
// Méthode basique : méthode de bissection 
// Défaut : Optimisation univariée, or plusieurs paramètres !!
// On fixe 3 paramètres sur 4 et on effectue la méthode sur le dernier paramètre
// Amélioration potentielle : Utiliser une méthode d'optimisation multivariée (descente de gradient?) 


double f(param_t* param, double param_test){
    // Fonction dont on cherche la racine
    double param_prev = param->l;
    param->l = param_test;
    printf("Testing for param : %.9e\n",param_test);
    double *f_actuals = compute_band_freq(0,param);
    // double *f_actuals = compute_band_freq(0,param);
    printf("Freq act = %f\n",f_actuals[0]);
    param->l = param_prev;
    return fabs(f_actuals[0] - f_target);
}

double bissection_method(double low_bound, double up_bound,param_t *param){
    printf("\n------------------------------------------------------\n");
    printf(" Bissection method !\n");
    printf("------------------------------------------------------\n");

    // printf("Checking if bounds have opposite signs\n");
    // if (f(param,low_bound) * f(param,up_bound) > 0){
    //     printf("f(low) and f(up must have opposite signs.");
    //     return -1; 
    // }
    
    double mid = (low_bound + up_bound) / 2.0;
    
    double err = fabs(f(param,mid));
    printf("--- > Initial error : %.9e\n",err);
    int iter =0;
    // while (fabs(f(param,mid)) > TOL && iter < 1000) {
     while (fabs(up_bound - low_bound) > TOL && iter < 1000) {
      if (f(param,mid) < TOL) {
        double err = fabs(f(param,mid)-f_target);
        printf("Freq computed ! Err final = %.9e\n!!",err);
        return mid;                     // Racine trouvée !
      }
      else if (f(param,mid)*f(param,low_bound) < 0 ){
        printf("Changing upper bound from %.9e to %.9e\n",up_bound,mid);
        up_bound = mid;
      }  
      else{
        printf("Changing lower bound from %.9e to %.9e\n",low_bound,mid);
        low_bound = mid;
      }
      mid = (up_bound + low_bound) / 2.0;
      printf("OK\n");
      // printf("Err act = %.9e\n",fabs(f(param,mid)));
      printf("low bound = %.9e\n",low_bound);
      printf("up bound = %.9e\n",up_bound);
      iter++;
    }
    return mid;
}

double f_half(param_t* param, double param_test){
    // Fonction dont on cherche la racine
    double param_prev = param->l;
    param->l = param_test;
    // printf("Testing for param : %.9e\n",param_test);
    double *f_actuals = compute_band_freq_half(0,param);
    // printf("Freq act = %f\n",f_actuals[0]);
    param->l = param_prev;
    return fabs(f_actuals[0] - f_target);
}

double bissection_half_method(double low_bound, double up_bound,param_t *param){
    printf("\n------------------------------------------------------\n");
    printf(" Bissection method !\n");
    printf("------------------------------------------------------\n");

    // printf("Checking if bounds have opposite signs\n");
    // if (f(param,low_bound) * f(param,up_bound) > 0){
    //     printf("f(low) and f(up must have opposite signs.");
    //     return -1; 
    // }
    
    double mid = (low_bound + up_bound) / 2.0;
    
    double err = fabs(f_half(param,mid));
    printf("--- > Initial error : %.9e\n",err);
    int iter =0;
    //while (err > TOL && iter < 1000) {
     while (fabs(up_bound - low_bound) > TOL && iter < 1000) {
      if (f_half(param,mid) < TOL) {
        double err = fabs(f_half(param,mid)-f_target);
        printf("Freq computed ! Err final = %.9e\n!!",err);
        return mid;                     // Racine trouvée !
      }
      else if (f_half(param,mid)*f_half(param,low_bound) < 0 ){
        printf("Changing upper bound from %.9e to %.9e\n",up_bound,mid);
        up_bound = mid;
      }  
      else{
        printf("Changing lower bound from %.9e to %.9e\n",low_bound,mid);
        low_bound = mid;
      }
      mid = (up_bound + low_bound) / 2.0;
      printf("OK\n");
      // printf("Err act = %.9e\n",fabs(f(param,mid)));
      printf("low bound = %.9e\n",low_bound);
      printf("up bound = %.9e\n",up_bound);
      iter++;
    }
    return mid;
}