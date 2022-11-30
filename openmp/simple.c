#include <stdio.h>

#define TIMES 10000

void compute(float a, float *b);
void add(float a, float *b);

int main()
{
	float a,b;
	int j;

	for (j=0;j<TIMES;j++){
		a=a+1;
		compute(a,&b);}

	printf("Result:%E\n", b);
}

void compute (float a, float *b)
{
	int i;

	for(i=0;i<TIMES;i++)
		add(a,b);
}

void add(float a, float *b)
{
	*b=a+(*b);
}

