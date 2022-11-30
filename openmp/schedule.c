#include <stdio.h>
#include <math.h>

int main()
{
    float A[400][400], tmp;
    int i,j,k;
    
#pragma omp parallel shared(A), private(i,j,k,tmp)
    {
#pragma omp for
        for(i=0;i<400;i++)
            for(j=0;j<400;j++)
                A[i][j]=i*j;
        
#pragma omp for
        for(i=0;i<400;i++)
            for(j=0;j<i;j++)
                for(k=0;k<1000;k++)
                {
                    tmp=sqrt(pow(A[i][j],j));
                    A[i][j]=pow(A[j][i],j)/i+tmp;
                }
    }
    
    printf("A[1][1]=%f,A[2][2]=%f,A[399][399]=%f\n",A[1][1],A[2][2],A[399][399]);
}



