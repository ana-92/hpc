#include <stdio.h>
#include <math.h>

int main()
{
    int i;
    double x, scale;
  
    scale=2.78;
    x=0.0;
#pragma omp parallel for private(x)
    for (i=0;i<1000;i++)
        x=x+sqrt(i*scale);


  printf("x=%f\n", x);

}

