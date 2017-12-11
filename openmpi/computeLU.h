#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>

void print_matrix(double *A, int rows, int columns);
int luDecomposeMPI(char *filename, double *l, double *u, int *dim);
