#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "../headers/matrix.h"

Matrix * allocate_matrix(int m, int n) {
	Matrix * mat = (Matrix*) malloc(sizeof(Matrix));
	mat->m = m, mat->n = n;
	mat->data = (double*) calloc(m*n,sizeof(double));
	mat->a = (double**) calloc(m,sizeof(double*));
	for(int i = 0; i < m; i++)
		mat->a[i] = mat->data+i*n;
	return mat;
}
BandMatrix * allocate_band_matrix(int m, int k){
    BandMatrix * band= malloc(sizeof(Matrix));
    band->m=m;
    band->k=k;

    band->data=(double*)calloc(m * (2*k+1), sizeof(double));
    band->a=(double**) calloc(m,sizeof(double*));
    for(int i=0; i<m; i++){
        band->a[i]=band->data+  k + 2*k*i;
    }
    return band;
}
void free_matrix(Matrix * mat) {
	free(mat->a);
	free(mat->data);
	free(mat);

}

void free_band_matrix(BandMatrix * mat){
    free(mat->data);
    free(mat->a);
    free(mat);
}


void print_vector(double * v, int n) {
	for(int i = 0; i < n; i++)
		printf("%.3e ", v[i]);
	printf("\n");
}

void print_matrix(Matrix * A) {
	for(int i = 0; i < A->m; i++)
		print_vector(A->a[i], A->n);
}

int is_symmetric(Matrix * K) {
	int symmetric = 1;
	for(int i = 0; i < K->m; i++)
		for(int j = i+1; j < K->n; j++)
			if(fabs((K->a[i][j] - K->a[j][i]) / K->a[i][j]) > 1e-12) {
				printf("%d %d\n", i, j);
				printf("%lf %lf\n", K->a[i][j], K->a[j][i]);
				symmetric = 0;
			}
	return symmetric;
}

typedef struct {
	int i;		// index
	double x,y;	// coordinates
} Node;

// Comparateur
int cmp(const void * a, const void * b) {
	
	Node * na = (Node *) a;
	Node * nb = (Node *) b;
	if (na->x > nb->x) return 1;
	else{
		if(na->y > nb->y) return 1;
	}
	return -1;
}

int cmp_y(const void * a, const void * b) {
	
	Node * na = (Node *) a;
	Node * nb = (Node *) b;
	if (na->y > nb->y) return 1;
	else{
		if(na->x > nb->x) return 1;
	}
	return -1;
}

int compute_permutation(int * perm, double * coord, int n_nodes) {

	// // We assume perm is allocated but not initialized
	// for(int i = 0; i < n_nodes; i++) 
	// 	perm[i] = i;

	// qsort_r(perm, n_nodes, sizeof(int), coord, cmp);

	// Create Node structs
	Node * nodes = malloc(n_nodes * sizeof(Node));
	for(int i = 0; i < n_nodes; i++) {
		nodes[i].i = i;
		nodes[i].x = coord[2*i];
		nodes[i].y = coord[2*i+1];
	}

	// Sort nodes
	//By comparing x
	//qsort(nodes, n_nodes, sizeof(Node), cmp);
	//By comparing y
	qsort(nodes, n_nodes, sizeof(Node), cmp_y);

	// Fetch permutation (we assume perm is allocated)
	for(int i = 0; i < n_nodes; i++)
		perm[i] = nodes[i].i;

	return 0;
}

double** permute_matrix(double** matrix, int* perm, int size) {
     double** new_matrix = (double**) malloc(size * sizeof(double*));
    for (int i = 0; i < size; i++) {
        new_matrix[i] = (double*) malloc(size * sizeof(double));
        for (int j = 0; j < size; j++) {
            new_matrix[i][j] = matrix[perm[i]][perm[j]];
        }
    }
    return new_matrix;
}

double* permute_vect(double* v,int* perm, int n_nodes) {

	double* v_perm = (double*) malloc(sizeof(double)*n_nodes);

	for (int i = 0; i < n_nodes; i++)
	{
		v_perm[i] = v[perm[i]];
		
	}
	return v_perm;
}

double* permute_vect_inv(double* v,int* perm, int size) {

	double* v_perm = (double*) malloc(sizeof(double)*size);

	for (int i = 0; i < size; i++)
	{
		v_perm[perm[i]] = v[i];
	}
	return v_perm;
}

int bandwidth(double **A, int m) {
    int i, j, k, bw = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            if (A[i][j] != 0) {
                k = abs(i-j);
                if (k > bw) {
                    bw = k;
                }
            }
        }
    }
    return bw;
}