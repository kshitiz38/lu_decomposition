#include <LUDecomp.h>
// -------------------------------------- MAIN -------------------------------------------------
int main(int argc, char * argv[])
{
    int i, j, k, n=5, ITERS=1;
    double* l;
    double* u;
    double* O;
    double begin, end, total_time=0.0, sum;

    double* A = new double[n*n];
    for(i=0; i<n*n; i++){
        A[i] = (double)(rand()%10);
    }

    printf("Original Matrix\n");
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            printf("%f ", A[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");

    for(j=0; j<ITERS; j++){

        l = new double[n*n]();
        u = new double[n*n]();
        begin = omp_get_wtime();
        i = luDecomposeP(A, l, u, n);
        end = omp_get_wtime();
        total_time += end - begin;
        if(i==0){
           printf("Return: Success!\n");
        } else {
           printf("Return: Failure!\n");
        }
    }


    printf("L Matrix\n");
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            printf("%f ", l[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("U Matrix\n");
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            printf("%f ", u[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");

    O = new double[n*n]();
    for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			sum=0;
			for(k=0; k<n; k++)
			{
				sum = sum + getValueAtIJ(l, n, i+1, k+1) * getValueAtIJ(u, n, k+1, j+1);
			}
            setValueAtIJ(O, n, i+1, j+1, sum);
		}
	}

    printf("LU Matrix Product\n");
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            printf("%f ", O[i*n+j]);
        }
        printf("\n");
    }
    printf("\n==================================================\n");
    printf("Time Taken for %d Iterations: %f, with N = %d", ITERS, total_time, n);
    printf("\n==================================================\n\n");
}
