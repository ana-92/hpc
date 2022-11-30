/**
 * Parallle programming class
 * Intermediate MPI, Excercise 12
 * Cannon’s algorithm for multiplying two matrix with non-blocking routines.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void fill_matrix(double *A, double *B, double *C, double *D, int size)
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
            D[i * size + j] = 0.0;
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

void blocking_routines(double *temporal, int *dims, double *buffer, double *A, double *B, double *C, int totalSize, int n_dim, int right, int left, int down, int up, MPI_Status status, MPI_Comm world_comm)
{
    for (int i = 0; i < dims[0]; i++)
    {
        matrix_mult(A, B, C, n_dim);

        if (i == dims[0] - 1)
        {
            break;
        }

        MPI_Send(A, totalSize, MPI_DOUBLE, left, 1, world_comm);
        MPI_Recv(buffer, totalSize, MPI_DOUBLE, right, 1, world_comm, &status);
        temporal = buffer;
        buffer = A;
        A = temporal;
        MPI_Send(B, totalSize, MPI_DOUBLE, up, 2, world_comm);
        MPI_Recv(buffer, totalSize, MPI_DOUBLE, down, 2, world_comm, &status);
        temporal = buffer;
        buffer = B;
        B = temporal;
    }
}
void non_blocking_routines(double *temporal, int *dims, double *buffer, double *A, double *B, double *C, int totalSize, int n_dim, int right, int left, int down, int up, MPI_Status status, MPI_Comm world_comm)
{
    MPI_Request reqs[4];

    for (int i = 0; i < dims[0]; i++)
    {
        matrix_mult(A, B, C, n_dim);

        if (i == dims[0] - 1)
            break;

        MPI_Isend(A, totalSize, MPI_DOUBLE, left, 1, world_comm, &reqs[0]);
        MPI_Irecv(buffer, totalSize, MPI_DOUBLE, right, 1, world_comm, &reqs[2]);
        temporal = buffer;
        buffer = A;
        A = temporal;

        MPI_Isend(B, totalSize, MPI_DOUBLE, up, 2, world_comm, &reqs[1]);
        MPI_Irecv(buffer, totalSize, MPI_DOUBLE, down, 2, world_comm, &reqs[3]);
        temporal = buffer;
        buffer = B;
        B = temporal;

        for (int j = 0; j < 4; j++)
            MPI_Wait(&reqs[j], &status);
    }
}

int main(int argc, char *argv[])
{
    MPI_Comm world_comm;
    MPI_Status status;
    int myid, n_dim, size, left, right, up, down, N;
    double *A, *B, *C, *D, *buffer1, *buffer2, *temporal1, *temporal2;
    double start, end, start_non_blocking, end_non_blocking;

    N = 16;

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
    D = (double *)malloc(totalSize * sizeof(double));
    buffer1 = (double *)malloc(totalSize * sizeof(double));
    buffer2 = (double *)malloc(totalSize * sizeof(double));

    fill_matrix(A, B, C, D, n_dim);

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
    blocking_routines(temporal1, dims, buffer1, A, B, C, totalSize, n_dim, right, left, down, up, status, world_comm);
    MPI_Barrier(world_comm);
    end = MPI_Wtime();

    start_non_blocking = MPI_Wtime();
    non_blocking_routines(temporal2, dims, buffer2, A, B, D, totalSize, n_dim, right, left, down, up, status, world_comm);
    end_non_blocking = MPI_Wtime();

    if (myid == 0)
    {
        printf("Matrix C:\n");
        print_matrix(C, n_dim);
        printf("\n");
        printf("Calculation time blocking routines: %.4fs\n", end - start);
        printf("----------------------------------------\n");
        printf("Matrix D:\n");
        print_matrix(D, n_dim);
        printf("\n");
        printf("Calculation time non-blocking routines: %.4fs\n", end_non_blocking - start_non_blocking);
    }

    MPI_Comm_free(&world_comm);
    free(A);
    free(B);
    free(C);
    MPI_Finalize();

    return 0;
}
