 #include <math.h>
 #include <stdlib.h>
 #include <complex.h>
 #include <stdio.h>
 #include <time.h>

//#include "/usr/include/lapacke.h"

#define dim 2                         
 
 
  main()
  {
  int i, j;
  // die auszuwertende matrix A mit komplexen komponenten
  complex double A[dim][dim];
   
  double D_XX = 3.;
  double OM = 1.; 
  
  double delta_detuning = 0.010;
  double detuning_start = -1.00;
  double detuning;
  int N = 2;
     
  char FILE_NAME[2048+1024];
  snprintf(FILE_NAME,2048+1024,"SUCK%.1f_OMEGA_%.1f_MAX_DET%.2f.dat",D_XX,OM,N*delta_detuning); 

  FILE *f_dat;
  f_dat=fopen(FILE_NAME,"w");
  
  int n;
  for (n=0;n<N;n++)
  {
  detuning = n*delta_detuning+detuning_start;  
  detuning = 0.;
  /*Die zu diagonalisierende Matrix Definition*/
  /* dim=2*/ 
  A[0][0]=3. +I*0.;  A[0][1]=-2. +I*0.;   
  A[1][0]=4. +I*0.;  A[1][1]=-1. +I*0.;   
  
  /*// dim=4 
  A[0][0]=0. +I*0.;  A[0][1]=OM       +I*0.;   A[0][2]=0.       +I*0.; A[0][3]=0.      +I*0.;
  A[1][0]=OM +I*0.;  A[1][1]=0.5*D_XX +I*0.;   A[1][2]=0.       +I*0.; A[1][3]=OM      +I*0.;
  A[2][0]=0. +I*0.;  A[2][1]=0.       +I*0.;   A[2][2]=0.5*D_XX +I*0.; A[2][3]=0.      +I*0.;
  A[3][0]=0. +I*0.;  A[3][1]=OM       +I*0.;   A[3][2]=0.       +I*0.; A[3][3]=detuning+I*0.; 
  /*Definitionsende*/
  
   /*Ausgabe des Inputs*/
  printf("Hier die eingegebene Matrix:\n");
  for (i=0; i<dim;i++)
  { for (j=0; j<dim; j++) 
	 { printf("%f\t+I %f\t", creal(A[i][j]),cimag(A[i][j]));}
       printf("\n");} 	 
//  */
  /* Schreibe Matrix mit komplexen komponenten so um, dass lapack sie versteht A --> AT */
  // die matrix, die dann lapack uebergeben wird
  double AT[2*dim*dim];                 

  for (i=0; i<dim; i++)          
  {  for(j=0; j<dim; j++)
      { AT[2*(j+dim*i)]=creal(A[j][i]);
        AT[2*(j+dim*i)+1]=cimag(A[j][i]);}  }
  /*Aufrufen der Lapack-Routinen mit Compiler-Befehl cc -llapack */
   /*Parameter der Eigenwert- und Eigenvektorberechnung*/
  complex double w[dim], vl[dim][dim], vr[dim][dim], b[dim], WORK[2*dim], RWORK[2*dim];
  int n=dim;
  char jobvl='N';
  char jobvr='V';
  int lda=dim;
  int ldvl=dim;
  int ldvr=dim;
  int lwork=2*dim;
  int ok;

  zgeev_(&jobvl, &jobvr,&n, AT, &lda, w, vl, &ldvl, vr, &ldvr, WORK, &lwork, RWORK, &ok);

  /*Ausgabe der Eigenwerte und Eigenvektoren*/
 if (dim==4)  fprintf(f_dat,"%.10f \t %.10f \t %.10f \t %.10f \t %.10f \n",detuning,creal(w[0]),creal(w[1]),creal(w[2]),creal(w[3]));  
 if (dim==2)  fprintf(f_dat,"%.10f \t %.10f+I%.10f \t %.10f+I%.10f \n",detuning,creal(w[0]),cimag(w[0]),creal(w[1]),cimag(w[1]));   
  }

   return(0);
  }
