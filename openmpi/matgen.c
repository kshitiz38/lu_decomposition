#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


int main(){
    int n = 10;
    double val;
    int i, j, k;

    // double matrix[9] = {1.0,2.0,4.0,3.0,8.0,14.0,2.0,6.0,13.0};

    FILE *fp;
    fp=fopen("matrix", "w");
    for(i=0; i<n; i++){
        for(j=0; j<n-1; j++){
            val = (double)(rand()%1000 - 500);
            fprintf(fp, "%f ", val);
        }
        val = (double)(rand()%1000 - 500);
        fprintf(fp, "%f\n", val);
    }
}
