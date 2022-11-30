/**
 * Parallle programming class
 * Intermediate MPI, Excercise 10
 * Cannon’s algorithm for multiplying two matrix.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void fill_matrix(double *A, double *B, double *C, int size)
{
    int seed = 0;
    srand(seed);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            // Fill with random nums
            A[i * size + j] = 5 - (int)(10.0 * rand() / (RAND_MAX + 1.0));
            B[i * size + j] = 5 - (int)(10.0 * rand() / (RAND_MAX + 1.0));
            // Initialize C
            C[i * size + j] = 0.0;
        }
    }
}

// Matrix multiplication
void matrix_mult(double *A, double *B, double *C, int size)
{
    for (int i = 0; i < size; i++)
        for (int k = 0; k < size; k++)
            for (int j = 0; j < size; j++)
                C[i * size + j] += A[i * size + k] * B[k * size + j];
}

void print_matrix(double *matrix, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("matrix[%d][%d] = %lf ", i, j, matrix[i * size + j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    MPI_Comm world_comm;
    MPI_Status status;
    int myid, n_dim, size, left, right, up, down, N;
    double *A, *B, *C, *buffer, *temporal;
    double start, end;

    N=16;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create Dimensions
    int dims[2];
    int periods[2];
    dims[0] = dims[1] = sqrt(size);
    int sq_int = (int)dims[0];
    periods[0] = periods[1] = 1;

    if (sq_int * sq_int != size)
    {
        if (myid == 0)
            printf("The number of processors must be perfect square.\n");
        MPI_Finalize();
        return 0;
    }

    MPI_Dims_create(size, 2, dims);
    // Local matrix dimension
    n_dim = N / dims[0];
    int totalSize = n_dim * n_dim;

    // Memory allocation
    A = (double *)malloc(totalSize * sizeof(double));
    B = (double *)malloc(totalSize * sizeof(double));
    C = (double *)malloc(totalSize * sizeof(double));
    buffer = (double *)malloc(totalSize * sizeof(double));

    fill_matrix(A, B, C, n_dim);

    if (myid == 0)
    {   
        printf("N_dim:%d\n", n_dim);
        printf("----------------------------------------\n");
        printf("Matrix A:\n");
        print_matrix(A, n_dim);
        printf("----------------------------------------\n");
        printf("Matrix B:\n");
        print_matrix(B, n_dim);
        printf("----------------------------------------\n");
    }

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &world_comm);
    MPI_Cart_shift(world_comm, 0, 1, &left, &right);
    MPI_Cart_shift(world_comm, 1, 1, &up, &down);

    start = MPI_Wtime();

    for (int i = 0; i < dims[0]; i++)
    {
        matrix_mult(A, B, C, n_dim);

        if (i == dims[0] - 1){
            break;
        }

        //Tambien se puede utilizar MPI_Send y MPI_Recv por separado
        MPI_Sendrecv(A, totalSize, MPI_DOUBLE, left, 1, buffer, totalSize, MPI_DOUBLE, right, 1, world_comm, &status);
        temporal = buffer;
        buffer = A;
        A = temporal;
        MPI_Sendrecv(B, totalSize, MPI_DOUBLE, up, 2, buffer, totalSize, MPI_DOUBLE, down, 2, world_comm, &status);
        temporal = buffer;
        buffer = B;
        B = temporal;
    }

    MPI_Barrier(world_comm);

    end = MPI_Wtime();

    if (myid == 0)
    {
        printf("Matrix C:\n");
        print_matrix(C, n_dim);
        printf("----------------------------------------\n");
        printf("\n");
        printf("Calculation time: %.4fs\n", end - start);
    }


    MPI_Comm_free(&world_comm);
    free(A);
    free(B);
    free(C);
    MPI_Finalize();

    return 0;
}
