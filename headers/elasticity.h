#ifndef _FEM_MATRICES_H_
#define _FEM_MATRICES_H_

#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include "../../gmsh-sdk/include/gmshc.h"
// #include "../../gmsh-4.11.1-Windows64-sdk/include/gmshc.h"
#include <string.h>

int assemble_system(Matrix** K, Matrix** M, double ** coord, size_t** boundary_nodes, size_t* n_boundary_nodes, double E, double nu, double rho);
int assemble_system_half(Matrix** K, Matrix** M, double** coord, size_t** boundary_nodes, size_t* n_boundary_nodes, double E, double nu, double rho,size_t** sym_nodes, size_t* n_sym_nodes);
void visualize_in_gmsh(double* SOL, int n_nodes);
void get_boundary_nodes(size_t** bnd_nTags, size_t* bnd_nTags_n);
void get_boundary_nodes_half(size_t** bnd_nTags, size_t* bnd_nTags_n, size_t** sym_nodes, size_t* n_sym_nodes);
#endif
