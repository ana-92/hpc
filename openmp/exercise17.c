#include <omp.h>
#include <stdio.h>

int main()
{
    int a[20], b[20], c[20], i;
    int d[20], e[20];
    
    for(i=0;i<20;i++)
    {
        a[i]=0;
        b[i]=i;
        c[i]=2*i;
        d[i]=19-i;
    }
#pragma omp parallel private(i)
    {
#pragma omp for nowait
        for(i=0;i<20;i++)
            a[i]=b[i]+c[i]+100*i;
#pragma omp for
        for(i=0;i<20;i++)
            e[i]=a[d[i]];
    }
    for(i=0;i<20;i++)
        printf("%d, %d, %d\n", i, e[i], a[d[i]]);
}
