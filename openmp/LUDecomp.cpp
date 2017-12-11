#include "./LUDecomp.h"

double getValueAtIJ(double *M, int n, int i, int j){
    return M[ (i-1)*n + j-1 ];
}

void setValueAtIJ(double *M, int n, int i, int j, double val){
    M[ (i-1)*n + j-1] = val;
}

void interChangeRows(double *A, int n, int a, int b){
    int i;
    for(i=0; i<n; i++){
        double tmp = getValueAtIJ(A, n, a, i+1);
        setValueAtIJ(A, n, a, i+1, getValueAtIJ(A, n, b, i+1));
        setValueAtIJ(A, n, b, i+1, tmp);
    }
}

int largestRowInColumn(double *A, int n, int a){
    int i;
    double max = fabs(getValueAtIJ(A, n, a, a));
    int row = -1;
    for(i=a; i<n; i++){
        double current = fabs(getValueAtIJ(A, n, i+1, a));
        if (current>max){
            row = i+1;
            max = current;
        }
    }
    return row;
}

void partialPivot(double *A, int n){
    int i;
    for(i=0; i<n; i++){
        int row = largestRowInColumn(A, n, i+1);
        if (row!=-1){
            interChangeRows(A, n, i+1, row);
        }
    }
}

int luDecompose(double *A, double *l, double *u, int n) {

    // partialPivot(A, n);
    int i, j, k, m, flag=0;
    for(i=0; i<n; i++){
        setValueAtIJ(l, n, i+1, i+1, 1.0);
        double sum_u=0.0, sum_l=0.0;
        for (k=0; k<i; k++){
            sum_u += getValueAtIJ(u, n, k+1, i+1) * getValueAtIJ(l, n, i+1, k+1);
        }
        setValueAtIJ(u, n, i+1, i+1, (getValueAtIJ(A, n, i+1, i+1) - sum_u));

        for(j=0; j<2*(n-i-1); j++){
            if(j<n-i-1){
                sum_u=0.0;
                for (k=0; k<i; k++){
                    sum_u += getValueAtIJ(u, n, k+1, j+i+2) * getValueAtIJ(l, n, i+1, k+1);
                }
                setValueAtIJ(u, n, i+1, j+i+2, (getValueAtIJ(A, n, i+1, j+i+2) - sum_u));
            }
            else {
                m = j%(n-i-1);
                sum_l=0.0;
                for (k=0; k<i; k++){
                    sum_l += getValueAtIJ(u, n, k+1, i+1) * getValueAtIJ(l, n, m+i+2, k+1);
                }
                double divd = getValueAtIJ(u, n, i+1, i+1);
                if(divd==0.0) flag=-1;
                setValueAtIJ(l, n, m+i+2, i+1, (getValueAtIJ(A, n, m+i+2, i+1) - sum_l)/divd);
            }
        }
    }
    return flag;
}

int luDecomposeP(double *A, double *l, double *u, int n) {

    // partialPivot(A, n);
    int i, j, m, k=0, flag=0;
    for(i=0; i<n; i++){
        setValueAtIJ(l, n, i+1, i+1, 1.0);
        double sum_u=0.0, sum_l=0.0, val=0.0, divd;

        #pragma omp parallel for reduction(+:sum_u) firstprivate(val, n, i) private(k)
        for (k=0; k<i; k++){
            sum_u += getValueAtIJ(u, n, k+1, i+1) * getValueAtIJ(l, n, i+1, k+1);
        }
        setValueAtIJ(u, n, i+1, i+1, (getValueAtIJ(A, n, i+1, i+1) - sum_u));

        #pragma omp parallel for firstprivate(n, i) private(sum_u, sum_l, k, divd, m)
        for(j=0; j<2*(n-i-1); j++){
            if(j<n-i-1){
                // #pragma omp parallel for reduction(+:sum_u) firstprivate(val, A, l, u, n, i, j) private(k)
                sum_u = 0.0;
                for (k=0; k<i; k++){
                    sum_u += getValueAtIJ(u, n, k+1, j+i+2) * getValueAtIJ(l, n, i+1, k+1);
                }
                setValueAtIJ(u, n, i+1, j+i+2, (getValueAtIJ(A, n, i+1, j+i+2) - sum_u));
            }
            else {
                m = j%(n-i-1);
                // #pragma omp parallel for reduction(+:sum_l) firstprivate(val, A, l, u, n, i, m) private(k)
                sum_l = 0.0;
                for (k=0; k<i; k++){
                    sum_l += getValueAtIJ(u, n, k+1, i+1) * getValueAtIJ(l, n, m+i+2, k+1);
                }
                divd = getValueAtIJ(u, n, i+1, i+1);
                if(divd==0.0) flag=-1;
                setValueAtIJ(l, n, m+i+2, i+1, (getValueAtIJ(A, n, m+i+2, i+1) - sum_l)/divd);
            }
        }
    }
    return flag;
}
// ------------------------------------- PARALLEL ---------------------------------------------
// int main(int argc, char * argv[]){}
//
// int main(int argc, char * argv[])
// {
//     int i, j, k, n=10, ITERS=1;
//     double* l;
//     double* u;
//     double* O;
//     double begin, end, total_time=0.0, sum;
//
//     double* A = new double[n*n];
//     for(i=0; i<n*n; i++){
//         A[i] = (double)(rand()%100);
//     }
//
//     printf("Original Matrix\n");
//     for(i=0; i<n; i++){
//         for(j=0; j<n; j++){
//             printf("%f ", A[i*n+j]);
//         }
//         printf("\n");
//     }
//     printf("\n");
//
//     for(j=0; j<ITERS; j++){
//
//         l = new double[n*n]();
//         u = new double[n*n]();
//         begin = omp_get_wtime();
//         i = luDecomposeP(A, l, u, n);
//         end = omp_get_wtime();
//         total_time += end - begin;
//         if(i==0){
//            printf("Return: Success!\n");
//         } else {
//            printf("Return: Failure!\n");
//         }
//     }
//
//
//     printf("L Matrix\n");
//     for(i=0; i<n; i++){
//         for(j=0; j<n; j++){
//             printf("%f ", l[i*n+j]);
//         }
//         printf("\n");
//     }
//     printf("\n");
//     printf("U Matrix\n");
//     for(i=0; i<n; i++){
//         for(j=0; j<n; j++){
//             printf("%f ", u[i*n+j]);
//         }
//         printf("\n");
//     }
//     printf("\n");
//
//     O = new double[n*n]();
//     for(i=0; i<n; i++)
// 	{
// 		for(j=0; j<n; j++)
// 		{
// 			sum=0;
// 			for(k=0; k<n; k++)
// 			{
// 				sum = sum + getValueAtIJ(l, n, i+1, k+1) * getValueAtIJ(u, n, k+1, j+1);
// 			}
//             setValueAtIJ(O, n, i+1, j+1, sum);
// 		}
// 	}
//
//     printf("LU Matrix Product\n");
//     for(i=0; i<n; i++){
//         for(j=0; j<n; j++){
//             printf("%f ", O[i*n+j]);
//         }
//         printf("\n");
//     }
//     printf("\n==================================================\n");
//     printf("Time Taken for %d Iterations: %f, with N = %d", ITERS, total_time, n);
//     printf("\n==================================================\n\n");
// }
