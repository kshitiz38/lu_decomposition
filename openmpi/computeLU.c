#include "computeLU.h"

double rawStringToDouble(const char *input) {
    const size_t bytesInDouble = 8;

    union {
        double value;
        unsigned char bytes[bytesInDouble];
    } u;

    unsigned char *output = u.bytes;

    int i;
    for(i = 0; i < bytesInDouble; ++i) {
        sscanf(input, "%02hhX", output);

        input += 2;
        ++output;
    }

    return u.value;
}

void print_matrix(double *A, int rows, int columns){
    int i, j;
    for(i=0; i<rows; i++){
        for(j=0; j<columns; j++){
            printf("%f ", A[i*columns+j]);
        }
        printf("\n");
    }
}

int luDecomposeMPI(char *filename, double *l, double *u, int *dim){

    int n = *dim;
    int pid, nprocs;
    double start, end, start_decomp, end_decomp;

    MPI_File fh;
    MPI_Status status;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(pid==0) start = MPI_Wtime();
    // printf("Rank is %d out of %d processors\n", pid, nprocs);

    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    int my_rows = n / nprocs;
    if (n%nprocs>pid) my_rows+=1;
    // printf("Processor %d has %d lines\n", pid, my_rows);

    int i=pid, j, k;
    const size_t bytesInDouble = 8;
    int line_size = bytesInDouble*2*n;
    double *A_local = (double *)malloc(my_rows*n*sizeof(double));
    char *read_buffer = (char*) malloc(line_size);
    while(i<n){

        MPI_File_seek(fh, i*(line_size), MPI_SEEK_SET);
        MPI_File_read(fh, read_buffer, line_size, MPI_BYTE, &status);
        // printf("Processor %d line - %s\n", pid, read_buffer);

        for(j=0; j<n; j++)
        {
            char tmp[bytesInDouble*2];
            strncpy(tmp, read_buffer+j*bytesInDouble*2, bytesInDouble*2);
            double val = rawStringToDouble(tmp);
            // printf("%f\n", val);
            A_local[(i/nprocs)*n+j] = val;
        }
        i += nprocs;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // if(pid==0) printf("\n######### Begin LUDecompisiton #########\n");
    if(pid==nprocs-1) start_decomp = MPI_Wtime();

    double *current_row = (double *)malloc(n*sizeof(double));
    for(i=0; i<n-1; i++){
        if(i%nprocs==pid){
            for(j=0; j<n; j++) {
                current_row[j] = A_local[(i/nprocs)*n+j];
            }
        }
        MPI_Bcast(current_row, n, MPI_DOUBLE, i%nprocs, MPI_COMM_WORLD);

        for(j=i+1; j<n; j++){
            if(j%nprocs==pid){
                if(current_row[i]==0.0) return -1;
                double pivot = A_local[(j/nprocs)*n+i] / current_row[i];
                for(k=i+1; k<n; k++){
                    A_local[(j/nprocs)*n+k] = A_local[(j/nprocs)*n+k] - (pivot * current_row[k]);
                }
                A_local[(j/nprocs)*n+i] = pivot;
            }
        }
    }
    if(pid==nprocs-1) end_decomp = MPI_Wtime();
    // if(pid==0) printf("######### LUDecompisiton Finished #########\n");

    // print_matrix(A_local, my_rows, n);

    double *current_buf = (double *)malloc(n*sizeof(double));
    double *current_l = (double *)malloc(n*sizeof(double));
    double *current_u = (double *)malloc(n*sizeof(double));

    FILE *lmatrix;
    FILE *umatrix;
    if(pid==0){
        lmatrix = fopen("lmatrix", "w");
        umatrix = fopen("umatrix", "w");
    }

    for(i=0; i<n; i++){
        if(i%nprocs==pid){
            for(j=0; j<n; j++) {
                current_buf[j] = A_local[(i/nprocs)*n+j];
            }
        }

        MPI_Bcast(current_buf, n, MPI_DOUBLE, i%nprocs, MPI_COMM_WORLD);
        if(pid==0){
            for(j=0; j<n; j++){
                if(j<i) current_l[j]=current_buf[j];
                else if(j==i) current_l[j]=1.0;
                else current_l[j]=0.0;
            }

            for(j=0; j<n; j++){
                if(j<i) current_u[j]=0.0;
                else current_u[j] = current_buf[j];
            }

            for(j=0; j<n; j++){
                if(j!=n-1) fprintf(lmatrix, "%f ", current_l[j]);
                else fprintf(lmatrix,"%f\n", current_l[j]);

            }

            for(j=0; j<n; j++){
                if(j!=n-1) fprintf(umatrix, "%f ", current_u[j]);
                else fprintf(umatrix,"%f\n", current_u[j]);

            }
        }
        // MPI_Barrier(MPI_COMM_WORLD);
    }
    if(pid==0){
        end = MPI_Wtime();
        printf("Total Time Taken Including File I/O (Read and Write Ops): %fs\n", (end - start));
    }
    if(pid==nprocs-1){
        printf("Total Time Taken Excluding I/O: %fs\n", (end_decomp - start_decomp));
    }
    MPI_Finalize();
    return 0;
}

int main(){
    int dim;
    FILE *fp = fopen("dims", "rb");
    if(fp){
        fscanf(fp, "%d", &dim);
        // printf("Size on Input Matrix is: %d\n", dim);
    }
    fclose(fp);
    double *l = NULL;
    double *u = NULL;
    int result = luDecomposeMPI("pmatrix", l, u, &dim);
    return 0;
}
