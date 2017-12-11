#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

double getValueAtIJ(double *M, int n, int i, int j);
void setValueAtIJ(double *M, int n, int i, int j, double val);
void interChangeRows(double *A, int n, int a, int b);
int largestRowInColumn(double *A, int n, int a);
void partialPivot(double *A, int n);
int luDecompose(double *A, double *l, double *u, int n);
int luDecomposeP(double *A, double *l, double *u, int n);
