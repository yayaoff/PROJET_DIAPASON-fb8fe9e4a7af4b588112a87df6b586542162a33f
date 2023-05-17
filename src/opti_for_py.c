#include "../headers/opti.h"

// Méthode basique : méthode de bissection 
// Défaut : Optimisation univariée, or plusieurs paramètres !!
// On fixe 3 paramètres sur 4 et on effectue la méthode sur le dernier paramètre
// Amélioration potentielle : Utiliser une méthode d'optimisation multivariée (descente de gradient?) 

// tag_vis == 1 si on souhaite visualiser le résultat , 0 sinon
double *compute_freq(int tag_vis,param_t *param){
    printf("Computing freq for params : \n r1 = %.9e, r2 = %.9e, l = %.9e, e = %.9e, meshSize = %.9e\n",
    param->r1, param->r2, param->l, param->e, param->meshSizeFactor);
  
  // Define physical constants
  double E = 0.7e11;  // Young's modulus for Aluminum
  double nu = 0.3;    // Poisson coefficient
  double rho = 3000;  // Density of Aluminum
  
  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(0, NULL, 0, 0, &ierr);
  
  designTuningFork(param->r1, param->r2, param->e, param->l, param->meshSizeFactor, param->filename);
  
  // Number of vibration modes to find
  int k = param->k;
  double *frequencies = malloc(sizeof(double)*k);

  // Assemble the 2 matrices of the linear elasticity problem: 
  // M is the mass matrix
  // K is the stiffness matrix
  Matrix *K, *M;
  size_t* boundary_nodes;
  size_t n_boundary_nodes;
  double * coord;
  assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);

  // Remove lines from matrix that are boundary
  Matrix *K_new;
  Matrix *M_new;
  remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL);

  // Compute M := K^{-1} M 
  int n = K_new->n;
  lu(K_new);
  double Mj[n]; // column j
  for (int j=0; j<n; j++){ // column by column
    for (int i=0; i<n; i++) Mj[i] = M_new->a[i][j];
    solve(K_new, Mj);
    for (int i=0; i<n; i++){
      M_new->a[i][j] = Mj[i];
    }
  } // -> Now, in M_new we have stored K^-1M

  // Power iteration + deflation to find k largest eigenvalues
  Matrix * A = M_new; // alias for conciseness
  double * v = calloc(A->m,sizeof(double));
  double lambda, freq;
  FILE * file = fopen(param->argv[2], "w"); // open file to write frequencies
  for(int ki = 0; ki < k; ki++) {
    lambda = power_iteration(A, v);
    // printf("lambda %d = %f\n",ki,lambda);
    freq = 1./(2*M_PI*sqrt(lambda));

    frequencies[ki] = freq;

    fprintf(file, "%.9lf ", freq);

    printf("lambda = %.9e, f = %.3lf\n", lambda, freq);
    // Deflate matrix
    for(int i = 0; i < A->m; i++)
      for(int j = 0; j < A->n; j++)
        A->a[i][j] -= lambda * v[i] * v[j];

    // Put appropriate BC and plot
    double * vall = calloc(K->m, sizeof(double));
    int iv = 0, i_bnd = 0; 
    for(int i = 0; i < K->m/2; i++) {
      if(i_bnd < n_boundary_nodes && i == boundary_nodes[i_bnd]) {
        i_bnd++;
        continue;
      }
      vall[2*(i)]   = v[2*iv];
      vall[2*(i)+1] = v[2*iv+1];
      iv++;
    }
    if(tag_vis ==1) visualize_in_gmsh(vall, K->m/2);
  }
  // fclose(file);
  if(tag_vis ==1) gmshFltkRun(&ierr);

  free_matrix (K);
  free_matrix (M);
  free_matrix (K_new);
  free_matrix (M_new);
  free(boundary_nodes);
  return frequencies;

}

// tag_vis == 1 si on souhaite visualiser le résultat , 0 sinon
double *compute_band_freq(int tag_vis,param_t *param){
    //printf("Computing freq for params : \n r1 = %.9e, r2 = %.9e, l = %.9e, e = %.9e, meshSize = %.9e\n",
    //param->r1, param->r2, param->l, param->e, param->meshSizeFactor);
  
  // Define physical constants
  double E = 0.7e11;  // Young's modulus for Aluminum
  double nu = 0.3;    // Poisson coefficient
  double rho = 3000;  // Density of Aluminum

  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(0, NULL, 0, 0, &ierr);
  
  designTuningFork(param->r1, param->r2, param->e, param->l, param->meshSizeFactor, NULL);
  
  // Number of vibration modes to find
  int k = param->k;
  double *frequencies = malloc(sizeof(double)*k);

  // Assemble the 2 matrices of the linear elasticity problem: 
  // M is the mass matrix
  // K is the stiffness matrix
  Matrix *K, *M;
  size_t* boundary_nodes;
  size_t n_boundary_nodes;
  double * coord;
  assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);

  //for Band solver
  BandMatrix *K_new;
  BandMatrix *M_new;
  remove_bnd_lines_Band(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL);

  //Compute permutations
  int n_nodes=K_new->m/2;
  int ** Adj=AdjacencyMatrix(K_new->a,K_new->m);
  int *perm=(int *) malloc(sizeof(int)* 2*n_nodes);
  for(int i=0; i<n_nodes; i++) perm[i]=i;
  reverseCuthillMcKee(Adj, 2*n_nodes, perm);
  
  K_new->a=permute_matrix(K_new->a,perm,2*n_nodes);
  //M_new->a=permute_matrix(M_new->a,perm,2*n_nodes);
  K_new->k=bandwidth(K_new->a,K_new->m);
  printf("Matrix dimension=%d \n", K_new->m);
  printf("bandwidth=%d \n",K_new->k);

  // Compute M := K^{-1} M 
  int n = K_new->m;
  //lu(K_new);
  lu_band(K_new);
  double *Mj=(double*) malloc(sizeof(double)*n); // column j
  for (int j=0; j<n; j++){ // column by column
    for (int i=0; i<n; i++) Mj[i] = M_new->a[i][j];
    Mj=permute_vect(Mj,perm,2*n_nodes);
    //solve(K_new, Mj);
    solve_band(K_new,Mj);
    Mj=permute_vect_inv(Mj,perm,2*n_nodes);
    for (int i=0; i<n; i++){
      M_new->a[i][j] = Mj[i];
    }
  } // -> Now, in M_new we have stored K^-1M
  printf("K-1M OK!\n");

  // Power iteration + deflation to find k largest eigenvalues
  BandMatrix * A = M_new; // alias for conciseness
  double * v = malloc(A->m * sizeof(double));
  double lambda, freq;
  FILE * file = fopen(param->argv[2], "w"); // open file to write frequencies
  for(int ki = 0; ki < k; ki++) {
    //lambda = power_iteration(A, v);
    lambda = power_iteration_band(A, v);
    freq = 1./(2*M_PI*sqrt(lambda));
    
    fprintf(file, "%.9lf ", freq);

    printf("lambda = %.9e, f = %.3lf\n", lambda, freq);
    frequencies[ki] = freq;
    
    // Deflate matrix
    for(int i = 0; i < A->m; i++)
      for(int j = 0; j < A->m; j++)
        A->a[i][j] -= lambda * v[i] * v[j];

    // Put appropriate BC and plot
    double * vall = calloc(K->m, sizeof(double));
    int iv = 0, i_bnd = 0; 
    for(int i = 0; i < K->m/2; i++) {
      if(i_bnd < n_boundary_nodes && i == boundary_nodes[i_bnd]) {
        i_bnd++;
        continue;
      }
      vall[2*(i)]   = v[2*iv];
      vall[2*(i)+1] = v[2*iv+1];
      iv++;
    }
    if(tag_vis ==1) visualize_in_gmsh(vall, K->m/2);
  }

  fclose(file);
  if(tag_vis ==1) gmshFltkRun(&ierr);

  free_matrix (K);
  free_matrix (M);
  free_band_matrix (K_new);
  free_band_matrix(M_new);
  free(boundary_nodes);
  return frequencies;
}

int main(int argc, char *argv[]){

    if (argc < 6){
    printf("Usage: \n"
			"./opti <k> <out> <r1> <r2> <e> <l> <out>\n" 
			"---------------------------- \n\n"
      "- k\n "
      "- out\n "
			"- r1\n "
			"- r2\n "
      "- e\n "
      "- l\n "
      "\n");
		return -1;
  } 
  param_t *params = malloc(sizeof(param_t));
  params->k = atoi(argv[1]); params->argc = argc; params->argv = argv; 
  params->r1 = atof(argv[3]) ; params->r2 = atof(argv[4]); params->e= atof(argv[5]);params->l = atof(argv[6]);
  params->meshSizeFactor=0.3;params->filename=argv[2];
  
  double *freqs = compute_band_freq(0,params);

  return 0;
}