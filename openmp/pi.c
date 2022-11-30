#include <stdio.h>
#include <math.h>

int main()
{

  double w,x,sum,pi,f,a;
  int n,i;

  n=1000000000;

  w=1.0/n;
  sum=0.0;

  for(i=0;i<n;i++){
     x=w*(i-0.5);
     sum=sum+4.0/(1.0+x*x);
  }

  pi=w*sum;

  printf("pi=%f\n", pi);

}
