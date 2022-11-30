#include <stdio.h>
/* Program for multiplication of D=A*B and E= A*C */

#define NRA 600 		/* number of rows in matrix A */
#define NCA 600        	/* number of columns in matrix A */
#define NCB 600		/* number of columns in matrix B */
#define NCC 600		/* number of columns in matrix C */

main()
{
int    i, j, k;			/* misc */
double a[NRA][NCA], 		/* matrix A to be multiplied */
       b[NCA][NCB],      	/* matrix B to be multiplied */
       c[NCA][NCC],		/* matrix C to be multiplied */
       d[NRA][NCB],      	/* result matrix A*B */
       e[NRA][NCC];		/* result matrix A*C */


   /* Initialize A, B, C, D and E matrices */
   for (i=0; i< NRA; i++)
      for (j=0; j< NCA; j++)
         a[i][j]= 1.;
   for (i=0; i< NCA; i++)
      for (j=0; j< NCB; j++)
         b[i][j]= 1.;
   for (i=0; i< NCA; i++)
      for (j=0; j< NCC; j++)
         c[i][j]= 2.;

   for(i=0;i< NRA;i++)
      for(j=0;j< NCB;j++)
         d[i][j] = 0.0;
   for(i=0;i< NRA;i++)
      for(j=0;j< NCC;j++)
         e[i][j] = 0.0;



   /* Perform matrix multiply A.B */
   for(i=0;i< NRA;i++)
      for(j=0;j< NCB;j++)
         for(k=0;k< NCA;k++)
            d[i][j]+= a[i][k] * b[k][j];


   /* Perform matrix multiply A.C */
   for(i=0;i< NRA;i++)
      for(j=0;j< NCC;j++)
         for(k=0;k< NCA;k++)
            e[i][j]+= a[i][k] * c[k][j];


   printf("Done\n");

   printf("d[0][0]= %f\n ", d[0][0]);
   printf("d[NRA-1][NCB-1]= %f\n ", d[NRA-1][NCB-1]);
   printf("e[0][0]= %f\n ", e[0][0]);
   printf("e[NRA-1][NCC-1]= %f\n ", e[NRA-1][NCC-1]);
   printf ("\n");
}

