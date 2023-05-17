#include <stdio.h>
#include "matrix.h"

void remove_bnd_lines (Matrix *K, Matrix *M, size_t *bnd_nodes, size_t n_bnd_nodes, Matrix **K_new, Matrix **M_new, int* invperm);


void remove_bnd_lines_Band(Matrix *K, Matrix *M, size_t *bnd_nodes, size_t n_bnd_nodes, BandMatrix**K_new, BandMatrix **M_new, int* invperm);
void remove_bnd_lines_half (Matrix *K, Matrix *M, size_t *bnd_nodes, size_t n_bnd_nodes, BandMatrix**K_new, BandMatrix **M_new, int* invperm,size_t *sym, size_t n_sym);

void remove_coord(double * coord, double * res ,size_t * boundary_nodes, size_t n_boundary_nodes, int n);

double power_iteration_band(BandMatrix *A, double *v) ;

// inner product of 2 vectors
double inner(double* a, double* b, int n);

// normalize vector 
double normalize(double* a, int n);

double power_iteration(Matrix *A, double *v);