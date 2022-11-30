#include <stdio.h>
/* Program for multiplication of D=A*B and E= A*C */

#define NRA 800 		/* number of rows in matrix A */
#define NCA 800        	/* number of columns in matrix A */
#define NCB 800		/* number of columns in matrix B */

int main()
{
int    i, j, k;			/* misc */
double a[NRA][NCA], 		/* matrix A to be multiplied */
       b[NCA][NCB],      	/* matrix B to be multiplied */
       d[NRA][NCB];      	/* result matrix A*B */


   /* Initialize A, B, D*/
   for (i=0; i< NRA; i++)
      for (j=0; j< NCA; j++)
         a[i][j]= 1.;
   for (i=0; i< NCA; i++)
      for (j=0; j< NCB; j++)
         b[i][j]= 1.;

   for(i=0;i< NRA;i++)
      for(j=0;j< NCB;j++)
         d[i][j] = 0.0;



   /* Perform matrix multiply A.B */
   for(i=0;i< NRA;i++)
      for(j=0;j< NCB;j++)
         for(k=0;k< NCA;k++)
            d[i][j]+= a[i][k] * b[k][j];


   printf("Done\n");

   printf("d[0][0]= %f\n ", d[0][0]);
   printf("d[NRA-1][NCB-1]= %f\n ", d[NRA-1][NCB-1]);
   printf ("\n");
}

