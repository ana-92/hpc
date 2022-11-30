/**
 * Author : Ana Izaguirre
 * Parallle programming class 
 * Intermediate MPI, Excercise 4
 * 
 * Use MPI Type indexed to create a derived datatype that corresponds to the
    lower triangular part of the following 4x4 matrix:
                             0  1  2  3  
                             4  5  6  7  
                             8  9  10 11
                             12 13 14 15
    the lower triangular part is the elements 0, 4, 5, 8, 9, 10, 12, 13, 14, 15. Process 0 should create
    the matrix as a one-dimensional array, create the derived datatype, and send the lower triangular
    part with a single call to MPI Send. Process 1 should receive the lower triangular part with a
    single call to MPI Recv and then print the data it received.
*/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>


void generate_matrix (int *A,int n){
    for(int i=0; i< n; i++){
        for(int j=0;j<n; j++){
                A[n*i+j] = n*i+j;
            }
        }
}

void initialize_matrix(int * result,int n){
    for(int i=0; i< n; i++){
        for(int j=0;j<n; j++){
                result[n*i+j] = 0;
            }
        }
}

void print_matrix (int *A,int n){
    for(int i=0; i< n; i++){
        for(int j=0;j<n; j++){
                printf("A[%d][%d] = %d ", i,j,A[n*i+j]);
            }
            printf("\n");
        }
    printf("\n");
}

void main (int argc, char *argv[]) {
    int N = 4;
    int A[N*N], result[N*N], blocklength[N], disp[N];
    int myid, size;
    MPI_Status status;
    MPI_Datatype T_triangle;


    /* Initialize MPI */
    MPI_Init(&argc, &argv); 

    /* Get my rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); 

    /* Get the total number of processors */
    MPI_Comm_size(MPI_COMM_WORLD, &size); 



    //Initialize matrix result
    initialize_matrix(result, N);

    // elem. per block and displacements
    for(int i=0; i<N; i++) {
        blocklength[i] = i+1;
        disp[i] = N*i;
    }

    //Create derived data type
    MPI_Type_indexed(N, blocklength, disp, MPI_INTEGER, &T_triangle);
    MPI_Type_commit(&T_triangle);


    if (myid == 0) {
        //Create the matrix as one dimensional
        generate_matrix(A,N);
        //Send
        MPI_Send(&A, 1, T_triangle, 1, 0, MPI_COMM_WORLD);
        //Print
        printf("Matrix A:\n");
        print_matrix(A,N);
        MPI_Type_free(&T_triangle);
    }
    else{
        //Receive
        MPI_Recv(result, 1, T_triangle, 0, 0, MPI_COMM_WORLD, &status);
        
        //Print
        printf("Matrix result:\n");
        print_matrix(result,N);
    }

  
    MPI_Finalize(); /* Terminate MPI */
}
