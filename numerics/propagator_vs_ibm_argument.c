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

 #define TIME_DYNAMICS_TO_NS 0.0200
 #define T_STEPS 12000
 int N=T_STEPS;  
 double delta_t = ((double)TIME_DYNAMICS_TO_NS) * 1000. * 1000. /((double)T_STEPS); 

 // ############# ELECTRONIC PARAMETERS (start) ########################################################
 #define GAM_R (0.000000633 *0.) // (fs)^(-1) ... T1 time 750ps=750 000fs, GAM_R = 1/2T1 = 0.000000666
 // (fs)^(-1) ... T1 time 600ps=600 000fs, GAM_R = 1/2T1 = 0.000000833
 // (fs)^(-1) ... T1 time 537ps=537 000fs, GAM_R = 1/2T1 = 0.000000933
 #define AW_VC (0.5)
 #define AW_CC (0.0)
 #define AW_VV (1.0)
 #define AW_CV (0.0)
 
 #define OmL      ( 10.005015 * 0.  ) 
 #define DELTA    ( 0.0135*0. )
 #define t_NULL   ( delta_t * N*0.4 )
 #define WIDTH    ( 1.05*1.*1000.) // FULL WIDTH AT HALF MAXIMUM .. paper

 #define PHON_FAK (1.0)
 // ############# PHONONIC PARAMETERS (start) ##########################################################
 #define PHONON_COUPLING_OUTPUT
 #define TEMPERATUR (4.0001)//(2.01)//(1.7500100001) 
 #define BOSE(x) ( 1./( exp( (x) /(kB*TEMPERATUR)) - 1. ))
 #define N_b ( 3000 ) //(2000)
 #define q_0 0.0000001
 #define q_N (2.0000)
 double delta_q = ((q_N-q_0)/((double)N_b));

 #define hbar_omega_lo 0.0364   // in eV
 #define rez_epsilon_strich 0.0119347
 // 1/e_strich = 1/e_high_freq - 1/e_stat
    #define epsilon_high_freq 10.9
    #define epsilon_static 12.53
    #define rez_epsilon_null   18.096272   // e^2/eV nm 
    #define el_mass 5.6856800 // fs2 eV/nm2
    #define hbar_omega_h  (0.015) // in eV
    #define eff_mass_h (0.45) 
   // Valenzband Confinement a_v=3.1 nm 
    #define eff_mass_el (0.063)
    #define hbar_omega_el (0.03)// in eV
    // Leitungsband Confinement a_c=5.8 nm 
    #define c_q (0.00511) // nm/fs
    // cla 5110 m/s ... 5110 10^9 nm / (10^15 fs) = 5110 10^-6 nm/fs = 0.00511 
    #define rho (33514170.0) // eVfs^2 /nm^5
    // kg = ... nutze Energie=Kraft x Weg , [kg] = 6.24151*10^3 ev fs^2/nm^5  
    // eV =1.6021766208(98)×10^19 J 
    // rho 5370 kg/m^3 --> 5370*6241,51 eV fs^2/nm^5 = 33514170 eV fs^2/nm^5
    #define D_c (-14.6) // eV
    #define D_v (-4.8)  // eV
 
 // ############# PHONONIC PARAMETERS (end) ############################################################
// ############# RK4 PARAMETERS (end) ###################################################################
 
		 
// #################### Parameter-Uebernahme Ende ############################
//////////////////////////////////////////////////////////////////////////////
// ################### normalized transition coupling element phonon #########

double gvc (double q)
 { double faktor,argument_h,argument_el;
   double temp_gvc; 
   argument_h  = -1.* hbar * q * hbar * q * 0.25 / (hbar_omega_h * eff_mass_h   * el_mass );
   argument_el = -1.* hbar * q * hbar * q * 0.25 / (hbar_omega_el * eff_mass_el * el_mass );
   faktor   = sqrt (0.5 * hbar * q / ( rho * c_q) );
   temp_gvc =  faktor * ( D_c * exp( argument_el) - D_v * exp( argument_h) ) * PHON_FAK;
   return temp_gvc;  
 } 

// #########################  MAIN    ######################################### 

 int main()                                     
 {  int k,i,m;
    double t;

  char date[256];  time_t curtime;   struct tm *loctime;  curtime = time (NULL);  loctime = localtime (&curtime);
  strftime (date, 256, "%d_%m", loctime); 
  char time[1024]; strftime(time,1024,"%H_%M_%S",loctime);
  
// ###########################################################################
// ########################## Coupling Element _ Rotating Frame LA ###########
// ###########################################################################

complex double *la_coup     = calloc(2*N_b,sizeof(double));
        double *n_ph        = calloc(  N_b,sizeof(double));

char FILE_NAME[2048+1024];
#ifdef PHONON_COUPLING_OUTPUT
FILE * f_coup;
snprintf(FILE_NAME,2048+1024,"%s_%iK_PAR_LA_phon_COUPLING_w_Nb=%i_OmL=%.3f_PHON=%.2f.dat",date,(int)TEMPERATUR,N_b,OmL,PHON_FAK); 
f_coup=fopen(FILE_NAME,"w");
#endif

double la_polaron = 0.;
double float_q = 0.;

for (i=0;i<N_b;i++)
{
  float_q      =q_0 + delta_q * i;
  la_coup[i]   = gvc(float_q);
  la_polaron  += (delta_q/(2.*PI*PI))*(float_q/c_q) * la_coup[i] * conj(la_coup[i]);
  n_ph[i] = BOSE((hbar*c_q*float_q)); 
  #ifdef PHONON_COUPLING_OUTPUT
  fprintf(f_coup,"%.10f \t %.10f \t %.10f \n",float_q,cabs(la_coup[i]),n_ph[i]);
  #endif
}
printf("q_max=%.10f -- Polaronshift=%.10f -- Re[phi(0)]=%.10f -- -- Im[phi(0)]=%.10f \n",float_q,la_polaron,creal(phi),cimag(phi));
#ifdef PHONON_COUPLING_OUTPUT
fclose(f_coup);
#endif


printf("%s \n",time);

FILE * f_dat;
snprintf(FILE_NAME,2048+1024,"%s_%s_%iK_IBM_Benchmark_PROPAGATOR_w_Nb=%i_OmL=%.3f_PHON=%.2f.dat",date,time,(int)TEMPERATUR,N_b,OmL,PHON_FAK); 
f_dat=fopen(FILE_NAME,"w");

complex double IBM_arg = 0.;
complex double IBM_arg_NUM =0.;    
	
// ===================================================================================
      for (k=0; k<N; k++) // for increases the index after the loop is completed, k=0 in the first step
      {  
       t=delta_t*k;
	   IBM_arg     = 0. + I*0.;
       IBM_arg_NUM = 0. + I*0.;
  
         for(m=0;m<k;m++)
          {
          for (i=0; i<N_b; i++)
          {
           float_q  = q_0 + delta_q * i;   
                  
           phi_nm  = (2.*n_ph[i]+1.) * cos( c_q*float_q*delta_t*(k-m) ) - I * sin( c_q*float_q*delta_t*(k-m) ) ; 

           IBM_arg_NUM += delta_t*delta_t*delta_q*phi_nm*la_coup[i]*la_coup[i]*(float_q*float_q/(2.*PI*PI));

           IBM_arg +=  delta_q * (1./(2.*M_PI*M_PI)) * la_coup[i]  * conj(la_coup[i]) *   I*float_q*t/c_q; 
           IBM_arg += -delta_q * (1./(2.*M_PI*M_PI)) * la_coup[i]  * conj(la_coup[i]) * ( (2.*n_ph[i]+1.)*( 1.-cos(c_q*float_q*t) ) + I*sin(c_q*float_q*t) ) /(c_q*c_q);
		 
          }
		  }
           fprintf(f_dat,"%.10f  \t %.10f \t %.10f \n",t/1000.,cabs(IBM_arg),cabs(IBM_arg_NUM) );
           printf("%.10f  \t %.10f \t %.10f \n",t/1000.,cabs(IBM_arg),cabs(IBM_arg_NUM) );
	  }    
	     
  fclose(f_dat);    
  
  free(la_coup);
  free(n_ph);

  return 0;
}

