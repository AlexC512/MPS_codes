/*************************************************** 
 *   Copyright (C) 2008 by Alexander Carmele   *
 *   221b Baker St, Marylebone, London W1.     *
 ***************************************************/

#include "parameter_set.h"
// Zur Definition der Zustaende
#define N_l (4) // Zahl der Level
#define rho(a,b)      ( *(derivates     + (a)*N_l + (b)  ))
#define rho_out(a,b)  ( *(derivates_out + (a)*N_l + (b)  ))
#define rho_stst(a,b) ( *(derivates_stst+ (a)*N_l + (b)  ))

#define Zahl_der_DGL      ( N_l*N_l + 10  )
int N_DGL=Zahl_der_DGL;

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
int a,b,k,n;
double t;

time_t curtime; struct tm *loctime; curtime = time(NULL); loctime = localtime (&curtime); 
int hour = loctime -> tm_hour;int minute = loctime -> tm_min;int second = loctime -> tm_sec;
int year = 1900+loctime -> tm_year;int month = loctime -> tm_mon;int day  = loctime -> tm_mday;
printf("Date: %d.%d.%.d -- Time: %d:%d:%d  \n",day,month+1,year,hour,minute,second);  

char FILE_NAME[2048+1024];
snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_SPECTRUM_MOLLOW_BARE_OmL_H%.0f_Delta_H%.0f_GamR_H%.0f_GamP_H%.0f.dat",year,month+1,day,hour,minute,second,OmL/GAM_R,DELTA/GAM_R,GAM_R,GAM_P/GAM_R); 

FILE *f_spectrum;
f_spectrum=fopen(FILE_NAME,"w");
fprintf(f_spectrum,"# Frequency of Signal [Mc] - Signal_incoherent - Signal_complete \n");

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
// ########################## TAU_DOMAIN SOLUTION ############################
// ########################## TAU_INITIAL_CONDITIONS #########################
complex double obs_inf = OUTPUT[5]+I*OUTPUT[6];
// dieser steady state wert wird vom spektrumsignal abgezogen, um den kohaerenten anteil
// abzuziehen
double obs_coh = creal(obs_inf*conj(obs_inf));
// save the t density matrix to create the new tau density matrix
for (n=0;n<N_DGL;n++) { derivates_stst[n] = derivates[n]; derivates[n]=0.;}

// INITIAL VALUES FOR THE DIPOLE SPECTRUM
OUTPUT[5] = OUTPUT[1]; // Excited State Density
OUTPUT[6] = 0.; 
for(a=0;a<N_l;a++){for(b=0;b<N_l;b++){
// < |1><0|(t) |a><b| (t+tau)>  
if (a==3)  rho(a,b) = 0.;
// -----------
if (a==2)  rho(a,b) = 0.;
// -----------    
if (a==1)  rho(a,b) = 0.;
// ---------  
if ((a==0)&&(b==3)) rho(a,b) = rho_stst(1,3);
if ((a==0)&&(b==2)) rho(a,b) = rho_stst(1,2);
if ((a==0)&&(b==1)) rho(a,b) = rho_stst(1,1);
if ((a==0)&&(b==0)) rho(a,b) = rho_stst(1,0);	
}}

complex double *correlation     = calloc(2*MAX_STEPS,sizeof(double));
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
 correlation[k]     = (OUTPUT[5]+I*OUTPUT[6] -1.*obs_coh);
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
 // ############ END OF 4 VECTOR RUNGE KUTTA #######################
 STEADY_STATE_TEST(derivates,derivates_stst,OUTPUT,Zahl_der_DGL);
   if (OUTPUT[99]<progress_per_stst)
   {
    progress_per_stst = progress_per_stst/10.;  
    printf("PARAMETER=%.10f - ",parameter);   
    printf("Re[SPEC_OP]=%.10f - ",creal(correlation[k]));   
    printf("Im[SPEC_OP]=%.10f - ",cimag(correlation[k]));
    printf("TRACE=%.10f - "      ,OUTPUT[0]);   
    printf("StSt=%.10f \n"       ,OUTPUT[99]);   
    }
   }  
 else {printf("ACHTUNG: STEADY NICHT ERREICHT!\n");OUTPUT[99]=-0.10000000001;}
} // end of tau integration
int k_MAX=k; // End of correlation array
double spectrum_inc[Nk]; double spectrum_sum[Nk]; // spectrum vector
// calculate and modify spectrum // substract offset in procedure
calculate_spectrum(correlation,spectrum_inc,spectrum_sum,k_MAX,obs_coh);
// calculate saturation of peaks around laser and polariton
// write into file
double plot_x;
for(k=0;k<Nk;k++) 
{
//plot_x = hbar*(-0.5*k_interval+k_middle+delta_k*k);//+1328.00;
//plot_x = (-0.5*k_interval+k_middle+delta_k*k)/(2.*OmL);//+1328.00;
plot_x = (-0.5*k_interval+k_middle+delta_k*k);//+1328.00;    
//if ( (plot_x>1325.00) && (plot_x<1335.00))
fprintf(f_spectrum,"%.10f \t %.10f \t %.10f \t %.10f \n",plot_x,spectrum_inc[k],parameter,parameter*parameter);
}
fprintf(f_spectrum,"\n");

// free given arrays
free(temp);free(derivates);free(derivates_stst);free(k1);free(k2);free(OUTPUT);free(correlation);

fclose(f_spectrum);
return 0;
}

