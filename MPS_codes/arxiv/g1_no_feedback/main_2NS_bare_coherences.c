/*************************************************** 
 *   Copyright (C) 2008 by Alexander Carmele   *
 *   221b Baker St, Marylebone, London W1.     *
 ***************************************************/
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <time.h>
#define hbar ( 0.658212196 )     // eVfs // meVps // 畫Vns
#define PI M_PI
// ################ ZEITRAUM #######################
#define TIME_DYNAMICS_TO ( 0.0010 )
#define T_STEPS 10000
#define TAU_FAKTOR (2)
int N     = T_STEPS;
int N_TAU = TAU_FAKTOR*T_STEPS;
int STEP_IS_SAVED= 10; // jeder ... Schritt wird in ein File geschrieben!!
int TIME_OUTPUT  =2000;
int TIME_OUTPUT_TAU  = 2000;
double delta_t = ((double)TIME_DYNAMICS_TO) * 1000. * 1000. /( (double)T_STEPS); 
// ----------------------------- Parameters in Units
#define GAM_R ( 0.002 * 1.0) // 1/ps -- 1/T1 = 500ps 
#define GAM_P ( 0.000 * 10. ) // 1/ps -- 1/T2* = 50ps 
// -----------------------------
//----------------------------------
#define OmL    ( 0.008 * 1. * 1.  ) 
// -------- Spectral properties ------------
#define DELTA    (   0.0 * OmL )
#define MAX_STEPS (1300000) // to reach steady state in tau
// -------------- parameter plot -----------------------------
// #########################################################
// ############# RK4 #######################################
// #########################################################

#define CALCULATE_TEMP(a,b,cb,dest,i) ({for(i=0;i<N_DGL;i++) (dest)[i]=(a)[i]+(cb)*(b)[i];})
#define ADD_VECTOR(dest,in,in_fak,i) ({for(i=0;i<N_DGL;i++) dest[i]=dest[i]+(in_fak)* (in)[i];})

#define SECHSTEL 0.1666666666666666666666667
#define DRITTEL  0.3333333333333333333333333
#define HALFTE   0.5

//#include "parameter_set.h"
// Zur Definition der Zustaende
#define N_l (4) // Zahl der Level
#define rho(a,b)      ( *(derivates     + (a)*N_l + (b)  ))
#define rho_out(a,b)  ( *(derivates_out + (a)*N_l + (b)  ))
#define rho_stst(a,b) ( *(derivates_stst+ (a)*N_l + (b)  ))

#define Zahl_der_DGL      ( N_l*N_l + 10  )
int N_DGL=Zahl_der_DGL;



// ###############################################################################
// ###############################################################################
void STEADY_STATE_TEST
(complex double *derivates,
 complex double *deriv_stst,
         double *OUTPUT,
	 int N
)
{
 int i;
 OUTPUT[99]=0.;
 for (i=0;i<N;i++)
 {
   OUTPUT[99] += cabs(derivates[i]-deriv_stst[i]);
   deriv_stst[i] = derivates[i];
 }
 return;
}
// ###############################################################################
// ###############################################################################

// #################################################################################
// ############## DGL system  ######################################################
// #################################################################################

void calculate_next_time_step 
(
complex double *derivates, 
complex double *derivates_out,  
        double t,
        double *OUTPUT
)
{   

  
int a,b;

for (a=0;a<100;a++) OUTPUT[a]=0.;


for (a=0;a<N_l;a++)
{
for (b=0;b<N_l;b++)
{
  rho_out(a,b) = 0. + I*0.;	
// ----------------------------------------------------------------------
// ------------- Left part - Density Matrix -----------------------------
// ----------------------------------------------------------------------
// ##### ground state couples via the laser and cavity mode to both exciton states   
if  (a==0)            {rho_out(a,b) += +I * OmL                         * rho(1,b) ;}
// ##### first exciton couples via the laser and cavity mode to the ground and biexciton state        
if  (a==1)            {rho_out(a,b) += -1.* GAM_R                      * rho(a,b)
                                       +I * DELTA                      * rho(a,b)
                                       +I * OmL                         * rho(0,b) ;}                                           
// ----------------------------------------------------------------------
// ------------- Right part - Density Matrix -----------------------------
// ----------------------------------------------------------------------
// ##### ground state couples via the laser and cavity mode to both exciton states   
if  (b==0)            {rho_out(a,b) += -I * OmL                     * rho(a,1) ;}
// ##### first exciton couples via the laser and cavity mode to the ground and biexciton state        
if  (b==1)            {rho_out(a,b) += -1.* GAM_R                  * rho(a,b)
                                       -I * DELTA                  * rho(a,b)
                                       -I * OmL                     * rho(a,0) ;}                                           
// ----------------------------------------------------------------------
// ------------- Lindblad Radiative Verlust  ----------------------------
// ----------------------------------------------------------------------
if ((a==0)&&(b==0)) rho_out(a,b) += +2.* GAM_R * rho(1,1);
// ----------------------------------------------------------------------
// ------------- Lindblad Pure dephasing  ----------------------------
// ----------------------------------------------------------------------
if ((a==1)&&(b==0)) rho_out(a,b) += -1.* GAM_P * rho(a,b);
if ((a==0)&&(b==1)) rho_out(a,b) += -1.* GAM_P * rho(a,b);
// ----------------------------------------------------------------------
// ------------- Expectation Values -------------------------------------
// ----------------------------------------------------------------------
if ( (a==b)              ) OUTPUT[0] += creal(rho(a,b));
if ( (a==b)  && (a==0)   ) OUTPUT[4] += creal(rho(a,b));
if ( (a==b)  && (a==1)   ) OUTPUT[1] += creal(rho(a,b));
if ( (a==b)  && (a==2)   ) OUTPUT[2] += creal(rho(a,b));
if ( (a==b)  && (a==3)   ) OUTPUT[3] += creal(rho(a,b));

if ( (a==0) && (b==1)   ) OUTPUT[5] += creal(rho(a,b));
if ( (a==0) && (b==1)   ) OUTPUT[6] += cimag(rho(a,b));

} //b
} //a

  return ;
}

// ###############################################################################
// ###############################################################################
// #########################  MAIN    ############################################ 
// ###############################################################################
// ###############################################################################
int main()                                     
{
int k,n;
double t;

time_t curtime; struct tm *loctime; curtime = time(NULL); loctime = localtime (&curtime); 
int hour = loctime -> tm_hour;int minute = loctime -> tm_min;int second = loctime -> tm_sec;
int year = 1900+loctime -> tm_year;int month = loctime -> tm_mon;int day  = loctime -> tm_mday;
printf("Date: %d.%d.%.d -- Time: %d:%d:%d  \n",day,month+1,year,hour,minute,second);  

double parameter = DELTA;

// ###############################################################################
// ############################## INITIAL CONDITIONS #############################
// ###############################################################################
complex double *derivates      = calloc(2*N_DGL,sizeof(double));
complex double *derivates_stst = calloc(2*N_DGL,sizeof(double));
// Um den Anfangswert zu speichern, mit dem die Zwischen-Anfangswerte berechnet werden
complex double *temp = calloc(2*N_DGL,sizeof(double));
// Die errechneten Zwischenfunktionswerte - die Steigungen
complex double *k1 = calloc(2*N_DGL,sizeof(double));
complex double *k2 = calloc(2*N_DGL,sizeof(double));
        double *OUTPUT = calloc(100,sizeof(double)); 
for (n=0;n<N_DGL;n++) derivates[n] = 0. + I* 0.;
for (n=0;n<100;n++)   OUTPUT[n]    =0.;
// set initial conditions: quantum in the ground state
rho(0,0) = 1. ; // QD in the Ground State with Zero Photons in the Cavity
OUTPUT[0] = 1.; // Trace=1
// ###########################################################################
// ########################## TIME _DOMAIN SOLUTION ##########################
// ###########################################################################
char FILE_NAME[2048+1024];
FILE *f_pop;
snprintf(FILE_NAME,2048+1024,"RHO_POPULATION_EE_OmL_H%.4f_Delta_H%.4f_GamR_H%.4f_GamP_H%.4f.dat",OmL/GAM_R,DELTA/GAM_R,GAM_R,GAM_P/GAM_R); 
f_pop=fopen(FILE_NAME,"w");

int i; 
k=0; double progress_per_stst = 0.1; OUTPUT[99] = 100.;
// begin t integration
while(OUTPUT[99]>0.0000000001)
{  
 k++;
 t=delta_t*k;
 // ############### 4 VECTOR RUNGE KUTTA ###########################
 calculate_next_time_step(derivates,k1,t            ,OUTPUT); CALCULATE_TEMP(derivates,k1,delta_t*HALFTE,temp,i);
 calculate_next_time_step(temp     ,k2,t+delta_t*0.5,OUTPUT); CALCULATE_TEMP(derivates,k2,delta_t*HALFTE,temp,i);
 ADD_VECTOR(k1,k2,2.0,i);
 calculate_next_time_step(temp     ,k2,t+delta_t*0.5,OUTPUT); CALCULATE_TEMP(derivates,k2,delta_t,temp,i);
 ADD_VECTOR(k1,k2,2.0,i);
 calculate_next_time_step(temp     ,k2,t+delta_t    ,OUTPUT); ADD_VECTOR(k1,k2,1.0,i); 
 ADD_VECTOR(derivates,k1,delta_t*SECHSTEL,i);
 // ############ END OF 4 VECTOR RUNGE KUTTA #######################
 STEADY_STATE_TEST(derivates,derivates_stst,OUTPUT,Zahl_der_DGL);
 if ( (k % STEP_IS_SAVED) == 0 ) fprintf(f_pop,"%.10f \t %.10f \n",k*delta_t*GAM_R*4.,creal(rho(1,1)));  
 if (OUTPUT[99]<progress_per_stst)
    {
     progress_per_stst = progress_per_stst/10.;  
     printf("PARAMETER=%.5f - ",parameter);   
     printf("HH=%.10f - "   ,OUTPUT[1]);
     printf("VV=%.10f - "   ,OUTPUT[2]);
     printf("BB=%.10f - "   ,OUTPUT[3]);
     printf("GG=%.10f - "   ,OUTPUT[4]);
     printf("TRACE=%.10f - "    ,OUTPUT[0]);
     printf("StSt_TEST=%.10f \n",OUTPUT[99]);
    }  
} // end of time integration
fclose(f_pop);
// ########################## TAU_DOMAIN SOLUTION ############################
// --- save the steady state density matrix for the different correlations
double gg = creal(rho(0,0));
double ee = creal(rho(1,1));
double ge = creal(rho(0,1));
double eg = creal(rho(1,0));
double norm;
// 
FILE *f_corr;
snprintf(FILE_NAME,2048+1024,"RHO_COHERENCE_EG_GE_OmL_H%.4f_Delta_H%.4f_GamR_H%.4f_GamP_H%.4f.dat",OmL/GAM_R,DELTA/GAM_R,GAM_R,GAM_P/GAM_R); 
f_corr=fopen(FILE_NAME,"w");
// ###########################################################################
// ########################## EGEG_GEGE CORRELATIONS #########################
// ###########################################################################
for (n=0;n<N_DGL;n++) { derivates[n]=0.;}
// --- normalization dependent on the correlation
norm = creal(ee*ee);
// save the t density matrix to create the new tau density matrix
// INITIAL VALUES FOR THE DIPOLE SPECTRUM
rho(0,1) = 0.;
rho(1,0) = ee;
rho(1,1) = 0.;
rho(0,0) = ge;

OUTPUT[99] = 100.; // steady state value
progress_per_stst = 0.1; // show progress every 1o percent
k=0; // set k back
// ########################## TAU_DOMAIN SOLUTION ############################
// integrate in tau until steady state 
while(OUTPUT[99]>0.0000000001)
{  
 k++;
 t=delta_t*k;
 // store all correlation points in array to Fourier Transform later
 // the Fourier Transform is not done here to substract a constant offset
 if (k<MAX_STEPS) 
 { 
 if ( (k % STEP_IS_SAVED) == 0 ) fprintf(f_corr,"%.10f \t %.10f \n",k*delta_t*GAM_R*4.,creal(rho(1,0)) );  
 // ############### 4 VECTOR RUNGE KUTTA ###########################
 calculate_next_time_step(derivates, k1 , t,OUTPUT); 
 CALCULATE_TEMP(derivates,k1,delta_t*HALFTE,temp,i);
 calculate_next_time_step(temp, k2, t+delta_t*0.5,OUTPUT); 
 CALCULATE_TEMP(derivates,k2,delta_t*HALFTE,temp,i);
 ADD_VECTOR(k1,k2,2.0,i);
 calculate_next_time_step(temp, k2, t+delta_t*0.5,OUTPUT); 
 CALCULATE_TEMP(derivates,k2,delta_t,temp,i);
 ADD_VECTOR(k1,k2,2.0,i);
 calculate_next_time_step(temp, k2, t+delta_t,OUTPUT);  
 ADD_VECTOR(k1,k2,1.0,i);
 ADD_VECTOR(derivates,k1,delta_t*SECHSTEL,i);
 // ############ eND OF 4 VECTOR RUNGE KUTTA #######################
 STEADY_STATE_TEST(derivates,derivates_stst,OUTPUT,Zahl_der_DGL);
   if (OUTPUT[99]<progress_per_stst)
   {
    progress_per_stst = progress_per_stst/10.;  
    printf("PARAMETER=%.10f - ",parameter);   
    printf("Correlation=%.10f - ",creal(rho(1,1))/norm);   
    printf("TRACE=%.10f - "      ,OUTPUT[0]);   
    printf("StSt=%.10f \n"       ,OUTPUT[99]);   
    }
   }  
 else {printf("ACHTUNG: STEADY NICHT ERREICHT!\n");OUTPUT[99]=-0.10000000001;}
} // end of tau integration
fclose(f_corr);



// free given arrays
free(temp);free(derivates);free(derivates_stst);free(k1);free(k2);free(OUTPUT);


return 0;
}

