/*************************************************** 
 *   Copyright (C) 2008 by Alexander Carmele   *
 *   221b Baker St, Marylebone, London W1.     *
 ***************************************************/

 #include <math.h>
 #include <stdlib.h>
 #include <complex.h>
 #include <stdio.h>
 #include <time.h>
 

 #define hbar 0.658212196      // eVfs 
 #define PI M_PI
 #define kB 0.0000861745 // eV/K

 #define TIME_DYNAMICS_TO_NS 10.00
 #define T_STEPS 30000
 int N=T_STEPS;  
 double delta_t = ((double)TIME_DYNAMICS_TO_NS) * 1000. * 1000. /((double)T_STEPS); 

 // certain: with increasing pump power decreasing wigner delay
 // but for a wide range of pump powers wigner delay remains the same
 // thus: pump power is not very well defined via the wigner signatures
 // scaling the dynamics with gamma_rad
 
 // ############# ELECTRONIC PARAMETERS (start) ########################################################
 #define GAM_R (0.000000833 *1.) // (fs)^(-1) ... T1 time 750ps=750 000fs, GAM_R = 1/2T1 = 0.000000666
 // (fs)^(-1) ... T1 time 600ps=600 000fs, GAM_R = 1/2T1 = 0.000000833
 // (fs)^(-1) ... T1 time 537ps=537 000fs, GAM_R = 1/2T1 = 0.000000933
 #define GAM_P (0.00000055000  *1.) // (fs)^(-1) ... 
 #define AW_RHOVC_SC 0.0
 #define AW_RHOCC_SC 0.0
 
 #define OmL    ( 0.20 * 1.  ) 
 #define t_NULL ( delta_t * N*0.4 )
 #define WIDTH  (1.05*1000.*1000.) // FULL WIDTH AT HALF MAXIMUM .. paper
// note, the pulse width is not optimal for the maximum wigner delay!!!
// however, the maximum wigner delay at this pulse width is still intensity independent!! 
// assuming a pulse width but not intensity dependent ... so for a fixed pulse width, we have 
// there is always a maximum wigner delay for a pulse width!!!!! and this is reduced by pure dephasing
// !!!!!!!!!!!! ATTENTION !!!! NO PHONONS !!! ATTENTION !!!
#define DELTA  ( 0.0000500 *1.0 ) // only important if no_phonons is switched on
#define NO_PHONONS 
#define PHON_FAK (0.00)
 // !!!!!!!!!!!! ATTENTION !!!! NO PHONONS !!! ATTENTION !!!
// ############# ELECTRONIC PARAMETERS (end) ##########################################################
#define Np ( 48 ) // es war 160
 // ############# PHONONIC PARAMETERS (start) ##########################################################
// #define PHONON_COUPLING_OUTPUT
 #define TEMPERATUR (1.3)//(2.01)//(1.7500100001) 
 #define BOSE(x) ( 1./( exp( (x) /(kB*TEMPERATUR)) - 1. ))
 #define N_b ( 2 ) //(2000)
 #define q_0 0.0000001
 #define q_N (2.0000)
 double delta_q = ((q_N-q_0)/((double)N_b));

 #define OPT_PHONON_SHIFT (0.0002225000 * PHON_FAK*1.0) // point of maximum delay (resonance)!!
 #define KAPPA_PH (0.0000005*1.) 

 #define hbar_omega_lo 0.0364   // in eV
 #define rez_epsilon_strich 0.0119347
 // 1/e_strich = 1/e_high_freq - 1/e_stat
 #define epsilon_high_freq 10.9
 #define epsilon_static 12.53
 #define rez_epsilon_null   18.096272   // e^2/eV nm 
 #define hbar_omega_el (0.04 ) //(0.03)//0.030 // 0.050   // in eV
 #define hbar_omega_h  (0.02 ) //(0.015) // (0.024) // in eV
 #define eff_mass_el (0.063)//(0.067)      
 #define eff_mass_h (0.45) // 0.45
 #define el_mass 5.6856800 // fs2 eV/nm2
 #define c_q (0.00511) // nm/fs
 #define rho (33508800.0) // eVfs^2 /nm^5
 #define D_c (-14.6*1.0*0.8) // eV
 #define D_v (-4.8*1.4 *0.8)  // eV
 // ############# PHONONIC PARAMETERS (end) ############################################################

// ############# RK4 PARAMETERS (start) ################################################################
#define VC          ( *(derivates       + 1 )) 
#define VC_OUT      ( *(derivates_out   + 1 )) 

#define CC          ( *(derivates       + 2 )) 
#define CC_OUT      ( *(derivates_out   + 2 )) 
 
#define VCB(m)      ( *(derivates       + 10 + (m) )) 
#define VCB_OUT(m)  ( *(derivates_out   + 10 + (m) )) 

#define VCBd(m)     ( *(derivates       + 11 + N_b + (m) )) 
#define VCBd_OUT(m) ( *(derivates_out   + 11 + N_b + (m) )) 

#define CCB(m)      ( *(derivates       + 12 + 2*N_b + (m) )) 
#define CCB_OUT(m)  ( *(derivates_out   + 12 + 2*N_b + (m) )) 

#define B(m)        ( *(derivates       + 13 + 3*N_b + (m) )) 
#define B_OUT(m)    ( *(derivates_out   + 13 + 3*N_b + (m) )) 
 
 
#define Zahl_der_DGL ( 20 + 4 * N_b  )
int N_DGL=Zahl_der_DGL;

#define CALCULATE_TEMP(a,b,cb,dest,i) ({for(i=0;i<N_DGL;i++) (dest)[i]=(a)[i]+(cb)*(b)[i];})
#define ADD_VECTOR(dest,in,in_fak,i) ({for(i=0;i<N_DGL;i++) dest[i]=dest[i]+(in_fak)* (in)[i];})
// ############# RK4 PARAMETERS (end) ###################################################################
 
	 
// #################### Parameter-Uebernahme Ende ############################
//////////////////////////////////////////////////////////////////////////////
// ################### normalized transition coupling element phonon #########

complex double gvc (double q)
 { double faktor,argument_h,argument_el;
   complex double temp_gvc; 
   argument_h  = -1.* hbar * q * hbar * q * 0.25 / (hbar_omega_h * eff_mass_h   * el_mass );
   argument_el = -1.* hbar * q * hbar * q * 0.25 / (hbar_omega_el * eff_mass_el * el_mass );
   faktor   = sqrt (0.5 * hbar * q / ( rho * c_q) );
   temp_gvc =  faktor * ( D_c * exp( argument_el) - D_v * exp( argument_h) ) * PHON_FAK;
   return temp_gvc;  
 } 

// ########################## DGL-System #####################################

void calculate_next_time_step 
(
 complex double* derivates, 
 complex double* derivates_out,  
         double  t, 
 complex double* coup_la, 
 complex double* rf_phonon,
         double  puls,
         double *n_ph,
         double detuning
)
{   

int m;
double         float_q  = q_0;

CC_OUT =  - 2. *  GAM_R        *   CC 
          - 2. * cimag (puls   *   VC          )
       ;

VC_OUT =  - 1. * (GAM_R+GAM_P) *   VC
          - I  * detuning      *   VC
          + I  * puls          * ( 2.*CC - 1.) 
       ;
#ifndef NO_PHONONS                
for (m=0;m<N_b;m++)
{
  float_q   = q_0 + delta_q * m;
  // phonon occupation 
  // hier zum Vergleich der Polarisationszerfall nach delta-puls anregung
  VC_OUT += -I*delta_q *float_q*float_q*      coup_la[m]  *  VCBd(m) *      rf_phonon[m] 
            -I*delta_q *float_q*float_q* conj(coup_la[m]) *  VCB(m)  * conj(rf_phonon[m])
          ;		
  // in 2nd order Born, Beitraege fuer jede Ordnung in (p)
  VCBd_OUT(m) = - 1.*(GAM_R    + GAM_P)           *          VCBd(m)
                - I * detuning                    *          VCBd(m)
                - 1.* KAPPA_PH                    *          VCBd(m)
                + I * puls                        * conj( 2.*CCB(m) - B(m) )  
                - I * n_ph[m]*conj(coup_la[m])    *     VC                    * conj(rf_phonon[m])
               ;

  VCB_OUT(m) = - 1.*(GAM_R    + GAM_P)          *      VCB(m)
               - I * detuning                   *      VCB(m)
               - 1.* KAPPA_PH                   *      VCB(m)
               + I * puls                       * ( 2.*CCB(m) - B(m) )        
               - I * (1.+ n_ph[m]) * coup_la[m] *      VC                     *      rf_phonon[m]  
             ;  
             
  CCB_OUT(m)  = -2.* GAM_R                      * CCB(m)
                -1.* KAPPA_PH                   * CCB(m)
                +I * puls                       * ( VCB(m) - conj( VCBd(m) ) )
                +I * n_ph[m] * conj(coup_la[m]) *   CC                         *      rf_phonon[m]  
              ;  
              
  B_OUT(m)    = -1.* KAPPA_PH                   * B(m)  
                -I * conj(coup_la[m])           * CC                           *      rf_phonon[m]  
              ;  
             
	}
#endif
	return ;
    }

// #########################  MAIN    ######################################### 

 int main()                                     
 {  int k,i;
    double t;
    double temp0,temp1,temp2;
    
    temp1 = 0.; temp2 = 0.;
// ###########################################################################
// ########################## Coupling Element _ Rotating Frame LA ###########
// ###########################################################################

complex double *la_coup     = calloc(2*N_b,sizeof(double));
        double *n_ph        = calloc(  N_b,sizeof(double));
complex double *rf_fak_1    = calloc(2*N_b,sizeof(double));
complex double *rf_fak_2    = calloc(2*N_b,sizeof(double));
complex double *rf_fak_dt   = calloc(2*N_b,sizeof(double));

char FILE_NAME[2048+1024];
#ifdef PHONON_COUPLING_OUTPUT
FILE * f_coup;
snprintf(FILE_NAME,2048+1024,"%iK_PAR_LA_phon_COUPLING_w_Nb=%i_OmL=%.3f_PHON=%.2f_GAM_PH=%.2f_GAM_RAD=%.10f.dat",(int)TEMPERATUR,N_b,OmL,PHON_FAK,GAM_P/GAM_R,GAM_R); 
f_coup=fopen(FILE_NAME,"w");
#endif

double la_polaron = 0.;
double float_q = 0.;

for (i=0;i<N_b;i++)
{
  float_q      =q_0 + delta_q * i;
  la_coup[i]   = gvc(float_q);
  la_polaron  += float_q * float_q * delta_q * la_coup[i] * conj(la_coup[i]) / (c_q * float_q);
  rf_fak_dt[i] = cexp( I * c_q * float_q * delta_t * 0.5 );
  rf_fak_1[i]  = 1. + I * 0.;
  rf_fak_2[i]  = 1. + I * 0.;
  n_ph[i] = BOSE((hbar*c_q*float_q)) ; 
  #ifdef PHONON_COUPLING_OUTPUT
  fprintf(f_coup,"%.10f \t %.10f \t %.10f \n",float_q,cabs(la_coup[i]),c_q*hbar*float_q);
  #endif
}

printf("q_max=%.10f -- Polaronshift=%.10f \n",float_q,la_polaron);
#ifdef PHONON_COUPLING_OUTPUT
fclose(f_coup);
#endif
// ###########################################################################
// ########################## TIME _DOMAIN SOLUTION ##########################
// ###########################################################################

int p;
double detuning = 0.0002225000; // point of maximum delay (resonance)!!
int delay_max = 0;
double wigner_delay_max =0.;
double wigner_delay[Np+10];
double CC_INT[Np+10];
double CC_INT_MAX=-10.;


for( p=0;p<Np;p++)
{   
 #ifdef NO_PHONONS
 detuning = DELTA*(-1.+2.*p/Np);
 #else   
 detuning = OPT_PHONON_SHIFT*(0.90+p*0.20/Np);
 #endif
 
    complex double *derivates = calloc(2*N_DGL,sizeof(double));
    // Die errechneten Zwischenfunktionswerte - die Steigungen
    complex double *k1 = calloc(2*N_DGL,sizeof(double));
    complex double *k2 = calloc(2*N_DGL,sizeof(double));
    // Um den Anfangswert zu speichern, mit dem die Zwischen-Anfangswerte berechnet werden
    complex double *temp = calloc(2*N_DGL,sizeof(double));
 
CC = AW_RHOCC_SC;
VC = 0. + AW_RHOVC_SC * I;

double puls;
double puls_max  =-10.;
double rhoee_max  =-10.;
double t_puls_max=-10.;
double t_rhoee_max=-10.;

CC_INT[p]=0.;

//Beginn der Integration von 0fs bis N*Schrittweite in fs                            
   for (k=0; k<N; k++)
     {  
       t=delta_t*k;
       temp0= t*0.001; // Wert in Picosekunden
       puls=(OmL/WIDTH)* exp(-1.0*(t-t_NULL)*(t-t_NULL)*1.386/(WIDTH*WIDTH));
       temp2   =creal(CC); 
       CC_INT[p] += creal(CC)*delta_t;
       
        if ( creal(CC)>rhoee_max) { rhoee_max = creal(CC); t_rhoee_max =t; }
        if ( puls>puls_max      ) { puls_max  = puls;      t_puls_max  =t; }

       // #####################################################
       // Uebergabe an die Runge-Kutta-Prozedur ###############
       // #####################################################

          calculate_next_time_step(derivates, k1 , t, la_coup, rf_fak_2,puls,n_ph,detuning); 
          CALCULATE_TEMP(derivates,k1,delta_t*0.5,temp,i);

  
          for (i=0; i<N_b; i++)
          {
           rf_fak_1[i] = rf_fak_2[i] * rf_fak_dt[i];
           rf_fak_2[i] = rf_fak_1[i] * rf_fak_dt[i]; 
          }

          calculate_next_time_step(temp, k2, t+delta_t*0.5, la_coup, rf_fak_1,puls,n_ph,detuning); 
 	      CALCULATE_TEMP(derivates,k2,delta_t*0.5,temp,i);
          ADD_VECTOR(k1,k2,2.0,i);
          
          calculate_next_time_step(temp, k2, t+delta_t*0.5, la_coup, rf_fak_1,puls,n_ph,detuning); 
	      CALCULATE_TEMP(derivates,k2,delta_t,temp,i);
          ADD_VECTOR(k1,k2,2.0,i);
          
          calculate_next_time_step(temp, k2, t+delta_t, la_coup, rf_fak_2,puls,n_ph,detuning); 
          ADD_VECTOR(k1,k2,1.0,i);
          ADD_VECTOR(derivates,k1,delta_t*0.166666667,i);
              
    }

 wigner_delay[p]=(t_rhoee_max-t_puls_max);   
 if ( wigner_delay[p]>wigner_delay_max ) {delay_max = p; wigner_delay_max=wigner_delay[p];}
 if ( CC_INT[p]>CC_INT_MAX) CC_INT_MAX=CC_INT[p];
 
 printf("%.2f Prozent. -- Detuning=%.10f -- Detuning_num=%.10f ",p*100./Np,(detuning-0.0002225000)*1000.*1000.*hbar,detuning); 
 printf("Normalized rhoee_max=%.10f \n",rhoee_max);
 printf("Wigner delay %.15f ps  \n",(t_rhoee_max-t_puls_max)*0.001);

  free(temp);
  free(derivates);
  free(k1);
  free(k2);
}
// ###################################################
// ######## begin write data #########################
// ###################################################
FILE *f_dat;
snprintf(FILE_NAME,2048+1024,"%.1fK_PAR_DETUNING_w_phon_Np=%i_Nb=%i_OmL=%.3f_PHON=%.2f_GAM_PH=%.2f_GAM_RAD=%.2f_Dc%.2f_Dv%.2f_confh%.4f_confe%.4f_phon_kap%.3f.dat",TEMPERATUR,Np,N_b,OmL,PHON_FAK,GAM_P*1000.*1000.,GAM_R*1000.*1000.,-1.*D_c,-1.*D_v,hbar_omega_h,hbar_omega_el,KAPPA_PH*1000.*1000.); 
f_dat=fopen(FILE_NAME,"w");

fprintf(f_dat,"# max delay detuning unshifted -- max_delay@ %i \n",delay_max);
for( p=0;p<Np;p++)
{
  #ifdef NO_PHONONS
  detuning = DELTA*(1.*p-1.*delay_max)/Np;
  #else   
  detuning = OPT_PHONON_SHIFT*(p-delay_max)*0.20/Np;
  #endif   
  fprintf(f_dat,"%.10f \t %.10f \t %.10f \n",detuning*1000.*1000.*hbar,wigner_delay[p]*0.001,wigner_delay_max*0.001*CC_INT[p]/CC_INT_MAX);
  }
fclose(f_dat);
// ###################################################
// ######## end write data ###########################
// ###################################################


// ###################################################
// ######## begin calculate resonant dynamics ########
// ###################################################
// now plot the time dynamics for the resonant case! where we have found the maximum wigner time delay

 #ifdef NO_PHONONS
 detuning = DELTA*(-1.+2.*delay_max/Np);
 #else   
 detuning = OPT_PHONON_SHIFT*(0.90+delay_max*0.20/Np);
 #endif
 
    complex double *derivates = calloc(2*N_DGL,sizeof(double));
    // Die errechneten Zwischenfunktionswerte - die Steigungen
    complex double *k1 = calloc(2*N_DGL,sizeof(double));
    complex double *k2 = calloc(2*N_DGL,sizeof(double));
    // Um den Anfangswert zu speichern, mit dem die Zwischen-Anfangswerte berechnet werden
    complex double *temp = calloc(2*N_DGL,sizeof(double));
 
CC = AW_RHOCC_SC;
VC = 0. + AW_RHOVC_SC * I;

double puls;
double puls_max  =-10.;
double rhoee_max  =-10.;
double t_puls_max=-10.;
double t_rhoee_max=-10.;

double CC_dyn[N];
complex double VC_dyn[N];
double PULS_dyn[N];


//Beginn der Integration von 0fs bis N*Schrittweite in fs                            
   for (k=0; k<N; k++)
     {  
       t=delta_t*k;
       temp0= t*0.001; // Wert in Picosekunden
       puls=(OmL/WIDTH)* exp(-1.0*(t-t_NULL)*(t-t_NULL)*1.386/(WIDTH*WIDTH));
       temp2   =creal(CC); 
       CC_dyn[k]   = creal(CC);
       VC_dyn[k]   = VC;
       PULS_dyn[k] = puls;
       
        if ( creal(CC)>rhoee_max) { rhoee_max = creal(CC); t_rhoee_max =t; }
        if ( puls>puls_max      ) { puls_max  = puls;      t_puls_max  =t; }

       // #####################################################
       // Uebergabe an die Runge-Kutta-Prozedur ###############
       // #####################################################

          calculate_next_time_step(derivates, k1 , t, la_coup, rf_fak_2,puls,n_ph,detuning); 
          CALCULATE_TEMP(derivates,k1,delta_t*0.5,temp,i);

  
          for (i=0; i<N_b; i++)
          {
           rf_fak_1[i] = rf_fak_2[i] * rf_fak_dt[i];
           rf_fak_2[i] = rf_fak_1[i] * rf_fak_dt[i]; 
          }

          calculate_next_time_step(temp, k2, t+delta_t*0.5, la_coup, rf_fak_1,puls,n_ph,detuning); 
 	      CALCULATE_TEMP(derivates,k2,delta_t*0.5,temp,i);
          ADD_VECTOR(k1,k2,2.0,i);
          
          calculate_next_time_step(temp, k2, t+delta_t*0.5, la_coup, rf_fak_1,puls,n_ph,detuning); 
	      CALCULATE_TEMP(derivates,k2,delta_t,temp,i);
          ADD_VECTOR(k1,k2,2.0,i);
          
          calculate_next_time_step(temp, k2, t+delta_t, la_coup, rf_fak_2,puls,n_ph,detuning); 
          ADD_VECTOR(k1,k2,1.0,i);
          ADD_VECTOR(derivates,k1,delta_t*0.166666667,i);
              
    }
FILE *f_dyn;
snprintf(FILE_NAME,2048+1024,"%.1fK_Dynamics_Detuning_%.8f_LA_phon_COUPLING_w_Nb=%i_OmL=%.3f_PHON=%.2f_GAM_PH=%.2f_GAM_RAD=%.10f.dat",TEMPERATUR,detuning,N_b,OmL,PHON_FAK,GAM_P/GAM_R,GAM_R); 
f_dyn=fopen(FILE_NAME,"w");
for (k=0; k<N; k++) 
fprintf(f_dyn,"%.10f \t %.10f \t %.10f \t %.10f \n",(k*delta_t-t_rhoee_max)*0.001*0.001,CC_dyn[k]/rhoee_max,PULS_dyn[k]/puls_max,creal(VC_dyn[k]*conj(VC_dyn[k]))/rhoee_max); 
fclose(f_dyn);

 wigner_delay[p]=(t_rhoee_max-t_puls_max);   
 if ( wigner_delay[p]>wigner_delay_max ) {delay_max = p; wigner_delay_max=wigner_delay[p];}
 if ( CC_INT[p]>CC_INT_MAX) CC_INT_MAX=CC_INT[p];
 
 printf("%.2f Prozent. -- Detuning=%.10f -- Detuning_num=%.10f ",p*100./Np,(detuning-0.0002225000)*1000.*1000.*hbar,detuning); 
 printf("Normalized rhoee_max=%.10f \n",rhoee_max);
 printf("Wigner delay %.15f ps  \n",(t_rhoee_max-t_puls_max)*0.001);

  free(temp);
  free(derivates);
  free(k1);
  free(k2);

  // ###################################################
  // ######## end calculate resonant dynamics ######## 
  // ###################################################

  
  free(rf_fak_1);
  free(rf_fak_2);
  free(rf_fak_dt);
  free(la_coup);
  free(n_ph);

  return 0;
}

