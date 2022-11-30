#include <stdio.h>
#include <math.h>
#include <omp.h>

#define MAX 10000
#define MAXITER 100000

double curra(int n)
{
    int i;
    double x = 0.0;

    // Algo para que curre un poco
    for(i=0;i<MAXITER;i++)
	x += sqrt((double) n)/(double) MAX;
 
   return x;
}

int main(int argc, char* argv[])
{
   int input=30, i=0, tid, nthreads;
   int vector[MAX];
   double result[MAX];

   for(i=0;i<MAX;i++)
	vector[i] = i+1;
   vector[MAX-1] = 0;


   i = 0;
   while(vector[i]){ // mientras haya curro...
   	   		    result[i] = curra(vector[i]);
				i++;}

   printf("Resultado para %d = %f \n ", vector[8], result[8]);
 
}


