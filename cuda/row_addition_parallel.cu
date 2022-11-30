
/* 
   Parámetros opcionales (en este orden): sumavectores  #size #blk
   #size: número de elementos en cada vector
   #blk: hilos por bloque CUDA
*/

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
const int SIZE = 5376;    // Número predeterm. de elementos en los vectores
const int CUDA_BLK = 32;  // Número predeterm. hilos por bloque

/* 
   Para medir el tiempo transcurrido (elapsed time):

   resnfo: tipo de dato definido para abstraer la métrica de recursos a usar
   timenfo: tipo de dato definido para abstraer la métrica de tiempo a usar

   timestamp: abstrae función usada para tomar las muestras del tiempo transcurrido

   printtime: abstrae función usada para imprimir el tiempo transcurrido

   void myElapsedtime(resnfo start, resnfo end, timenfo *t): función para obtener 
   el tiempo transcurrido entre dos medidas
*/

#include <sys/time.h>
#include <sys/resource.h>

#ifdef _noWALL_
typedef struct rusage resnfo;
typedef struct _timenfo {
  double time;
  double systime;
} timenfo;
#define timestamp(sample) getrusage(RUSAGE_SELF, (sample))
#define printtime(t) printf("%15f s (%f user + %f sys) ",		\
			    t.time + t.systime, t.time, t.systime);
#else
typedef struct timeval resnfo;
typedef double timenfo;
#define timestamp(sample)     gettimeofday((sample), 0)
#define printtime(t) printf("%15f s ", t);
#endif

void myElapsedtime(const resnfo start, const resnfo end, timenfo *const t)
{
#ifdef _noWALL_
  t->time = (end.ru_utime.tv_sec + (end.ru_utime.tv_usec * 1E-6)) 
    - (start.ru_utime.tv_sec + (start.ru_utime.tv_usec * 1E-6));
  t->systime = (end.ru_stime.tv_sec + (end.ru_stime.tv_usec * 1E-6)) 
    - (start.ru_stime.tv_sec + (start.ru_stime.tv_usec * 1E-6));
#else
  *t = (end.tv_sec + (end.tv_usec * 1E-6)) 
    - (start.tv_sec + (start.tv_usec * 1E-6));
#endif /*_noWALL_*/
}



/**
 * Generar un vector con valores flotantes aleatorios
*/
void generar_vector( float *A, int N){
    srand(time(NULL));
    //Fue necesario linealizar el arreglo a A[nxn]
    for(int i=0; i< N; i++){
            for(int j=0;j<N; j++){
                A[N*i+j] = (rand( ) % 7501 ) / 1000.0f;
            }
        }
}

/**
 * Imprimir el valor de los vectores resultantes
*/
void imprimir_vectores (float *A,float *sum_seq, float *sum, int N){
    printf("\n");
    printf("****** Vector Result CPU ******\n");	
	for (int i = 0; i < N; ++i)
		printf("sum[%d] = %lf\n", i, sum_seq[i]);

    printf("\n");
    printf("****** Vector Result GPU ******\n");	
	for (int i = 0; i < N; ++i)
		printf("sum[%d] = %lf\n", i, sum[i]);
    
    printf("\n");

}

/**
 * Verificar si el vector generado en el CPU es igual al del GPU
*/
void comparar_vectores(float *sum_seq, float *sum, int N){
    int contador = 0;
    for (int i = 0; i < N; ++i){
        contador++;
        if(sum[i]!= sum_seq[i]){
            break;
        }
    }

    
    if(contador == N){
         printf("Los valores en cada arreglo son iguales\n");
    }
    else{
         printf("Los valores en cada arreglo NO son iguales\n");
    }   

}

/**
 * La suma del vector en el CPU
*/

void rowSums_seq(float* A, float* sum, int N){
    for(int i=0; i< N; i++){
        for(int j=0;j<N; j++){
                sum[i] += A[N*i+j];
            }
            
    }
}

/*
  Definición de nuestro kernel en CUDA
*/
__global__ void rowSums(float* A, float* sum, int N)
{
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    float tmpSum = 0;

    if (row < N){
        for (int k = 0; k< N; k++) {
            tmpSum += A[row * N + k] ;
        }
    }

    sum[row] = tmpSum;
              
 } 

int main(int argc, char *argv[])
{
    float *cA ,*cSum;
    
    // Para medir tiempos
    resnfo start, end;
    timenfo time;

    // Si existe mas de un argumento lo usamos sino utilizamos el default.
    int size = (argc > 1)?atoi (argv[1]):SIZE;
    int blk_size = (argc > 2)?atoi (argv[2]):CUDA_BLK;
    
    printf("Elementos en los arreglos =%d \n",size);
    printf("Hilos por bloque =%d \n",blk_size);

    int numBytesA = size*size*sizeof(float);
    int numBytesSum = size*sizeof(float);


    //Reserva de memoria en el CPU
    float *A = (float *) malloc(numBytesA);
    float *sum = (float *) malloc(numBytesSum);
    float *sum_seq = (float *) malloc(numBytesSum);

    //Generar vector de tamaño size*size (CPU)
    generar_vector(A,size);

    //Realizar la suma en el CPU
    timestamp(&start);
    rowSums_seq(A,sum_seq,size);
    timestamp(&end);

    myElapsedtime(start, end, &time);
    printtime(time);
    printf(" -> Sumar vectores en CPU \n");


    //Reserva de memoria en el GPU
    cudaMalloc((void**)&cA, numBytesA);
    cudaMalloc((void **) &cSum , numBytesSum) ;

    // CPU -> GPU
    cudaMemcpy(cA, A, numBytesA, cudaMemcpyHostToDevice);

    //Inicializar arreglo de sumas
    cudaMemset(cSum, 0, numBytesSum);

    // Bloque unidimensional de hilos (*blk_size* hilos)
    dim3 dimBlock(blk_size);

    // Rejilla unidimensional (*ceil(n/blk_size)* bloques)
    dim3 dimGrid((size + dimBlock.x - 1) / dimBlock.x);


	//Cálculo de la suma del vector en el GPU
    timestamp(&start);
    rowSums<<<dimGrid, dimBlock >>>(cA, cSum, size);
    timestamp(&end);

    myElapsedtime(start, end, &time);
    printtime(time);
    printf(" -> Sumar vectores en GPU \n");

    // GPU -> CPU
    cudaMemcpy(sum, cSum, numBytesSum, cudaMemcpyDeviceToHost); 

     //Imprimir vectores
    //imprimir_vectores(A,sum_seq,sum,size);

    //Comparar vectores
    comparar_vectores(sum_seq,sum,size);

	cudaFree(cA);
	cudaFree(cSum);

	return 0;
}