#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmsh-4.11.1-Windows64-sdk/include/gmshc.h"
#include "headers/matrix.h"
#include "headers/elasticity.h"
#include "math.h"
#include "headers/lu.h"
#include "headers/design.h"
#include "headers/eigen.h"
#include "headers/RCM.h"

 #define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })


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
  // Define physical constants
  double E = 0.7e11;  // Young's modulus for Aluminum
  double nu = 0.3;    // Poisson coefficient
  double rho = 3000;  // Density of Aluminum

	// Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);

  designTuningFork(6e-3, 11e-3, 38e-3, 82e-3, 0.3, NULL);
  
  // Number of vibration modes to find
  int k = atoi(argv[1]);

  // Assemble the 2 matrices of the linear elasticity problem: 
  // M is the mass matrix
  // K is the stiffness matrix
  Matrix *K, *M;
  size_t* boundary_nodes;
  size_t n_boundary_nodes;
  size_t* sym_nodes;
  size_t n_sym_nodes;
  double * coord;
  assemble_system_half(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho,&sym_nodes,&n_sym_nodes);
  //n_boundary_nodes++;
  // size_t* bnd_nodes=malloc(n_boundary_nodes*sizeof(size_t));

  // for(int i=0; i<n_boundary_nodes; i++){
  //   // if(boundary_nodes[i]==9){
  //   //   //printf("OK!\n");
  //   //   bnd_nodes[i]=3;
  //   //   continue;
  //   // }
  //   bnd_nodes[i]=boundary_nodes[i];
  // }
  // //bnd_nodes[n_boundary_nodes-1]=35;

  // //n_sym_nodes--;
  // size_t* sym=malloc(n_sym_nodes*sizeof(size_t));
  // int offset=0;
  // for(int i=0; i<n_sym_nodes; i++){
  //   // if(sym_nodes[i]==8){
  //   //   sym[offset]=9;
  //   //   offset++;
  //   //   continue;
  //   // }
  //   //if(sym_nodes[i]==164)continue;
  //   sym[offset]=sym_nodes[i];
  //   offset++;
  // }
  // //sym[offset]=192;

  //for Band solver
  BandMatrix *K_new;
  BandMatrix *M_new;
  remove_bnd_lines_half(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL, sym_nodes,n_sym_nodes);
  //Compute permutations
  int n_nodes=K_new->m;
  int ** Adj=AdjacencyMatrix(K_new->a,K_new->m);
  int *perm=(int *) malloc(sizeof(int)*n_nodes);
  for(int i=0; i<n_nodes; i++) perm[i]=i;
  reverseCuthillMcKee(Adj, n_nodes, perm);
  
  K_new->a=permute_matrix(K_new->a,perm,n_nodes);
  //M_new->a=permute_matrix(M_new->a,perm,n_nodes);
  K_new->k=bandwidth(K_new->a,K_new->m);
  printf("Matrix dimension=%d \n", K_new->m);
  printf("bandwidth=%d \n",K_new->k);

  // Matrix *K_new;
  // Matrix *M_new;
  // remove_bnd_lines(K, M, bnd_nodes, n_boundary_nodes, &K_new, &M_new, NULL, sym,n_sym_nodes);
  

  // Compute M := K^{-1} M 
  int n = K_new->m;
  //lu(K_new);
  lu_band(K_new);
  double *Mj=(double*) malloc(sizeof(double)*n); // column j
  for (int j=0; j<n; j++){ // column by column
    for (int i=0; i<n; i++) Mj[i] = M_new->a[i][j];
    Mj=permute_vect(Mj,perm,n_nodes);
    //solve(K_new, Mj);
    solve_band(K_new,Mj);
    Mj=permute_vect_inv(Mj,perm,n_nodes);
    for (int i=0; i<n; i++){
      M_new->a[i][j] = Mj[i];
    }
  } // -> Now, in M_new we have stored K^-1M
  printf("K-1M OK!\n");

  // Power iteration + deflation to find k largest eigenvalues
  BandMatrix * A = M_new; // alias for conciseness
  //Matrix * A = M_new;
  double * v = malloc(A->m * sizeof(double));
  double lambda, freq;
  FILE * file = fopen(argv[2], "w"); // open file to write frequencies
  for(int ki = 0; ki < k; ki++) {
    //lambda = power_iteration(A, v);
    lambda = power_iteration_band(A, v);
    freq = 1./(2*M_PI*sqrt(lambda));

    fprintf(file, "%.9lf ", freq);

    printf("lambda = %.9e, f = %.3lf\n", lambda, freq);
    // Deflate matrix
    for(int i = 0; i < A->m; i++)
      for(int j = 0; j < A->m; j++)
        A->a[i][j] -= lambda * v[i] * v[j];

    // Put appropriate BC and plot
    double * vall = calloc(K->m, sizeof(double));
    int iv = 0, i_bnd = 0, i_sym = 0, offset = 0; 
    for(int i = 0; i < K->m/2; i++) {
      if(i_bnd < n_boundary_nodes && i == boundary_nodes[i_bnd]) {
        if(i==0) i_sym++;
        i_bnd++;
        continue;
      }
      if (i_sym < n_sym_nodes && i == sym_nodes[i_sym])
      {
        vall[2*i] = v[2*iv +offset ];
        
        i_sym++;
        offset++;
        
        continue;
      }
      vall[2*(i)]   = v[2*iv+offset ];
      vall[2*(i)+1] = v[2*iv+1+offset ];
      iv++;
    }
    visualize_in_gmsh(vall, K->m/2);
  }

  fclose(file);

  gmshFltkRun(&ierr);

  free_matrix (K);
  free_matrix (M);
  free_band_matrix (K_new);
  free_band_matrix(M_new);
  // free_matrix (K_new);
  // free_matrix(M_new);
  free(boundary_nodes);

  return 0;
}
