#include <stdio.h>

int main()
{
  int i,j;
  double A[4][50];

  for (i=0;i<4;i++)
#pragma omp parallel for shared(A), private(i,j)
    for(j=0;j<50;j++)
      A[i][j]=i*j;

  printf("A[1][1]=%f,A[2][2]=%f,A[3][3]=%f\n", A[1][1],A[2][2],A[3][3]);

}

