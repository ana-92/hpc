#include <omp.h>
#include <math.h>
#include <stdio.h>

#define N 10000

int main()
{
    int index[N], it, i, work[N];
    float x[N];
    

    for(i=0;i<N;i++)
    {
    
        work[i]=i;
        index[i]=i%N;
        x[i]=(float)i;
    }

    for (it=0;it<100;it++)
    {
#pragma omp parallel for private(i)
        for(i=0;i<N;i++){
#pragma omp critical
            x[index[i]]=x[index[i]]+sqrt(pow(work[i],3));
    }
    }
    printf("%f, %f\n", x[6],x[100]);
}
