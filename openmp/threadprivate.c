#include <stdio.h>

int j;
#pragma omp threadprivate(j)

int main()
{
    j=1;
    
#pragma omp parallel
{
#pragma omp master
    {
        j=2;
    }
}
printf("j=%d\n",j);
}
