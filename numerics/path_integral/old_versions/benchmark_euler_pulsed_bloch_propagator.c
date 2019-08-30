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

 #define TIME_DYNAMICS_TO_NS 0.05000
 #define T_STEPS 150
 #define SHOW_EVERY_STEP (T_STEPS/10)
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
 
 #define OmL      ( 0.0000784 * 2.0  ) 
 #define DELTA    ( 0.0135*0. )
 #define t_NULL   ( delta_t * N*0.33 )
 #define WIDTH    ( 1.05*1.*22000.) // FULL WIDTH AT HALF MAXIMUM .. paper

 #define PHON_FAK (2.0)
 // ############# PHONONIC PARAMETERS (start) ##########################################################
// #define PHONON_COUPLING_OUTPUT
 #define TEMPERATUR (50.0001)//(2.01)//(1.7500100001) 
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


#define M_Om(m,n)    ( *(Hsys      + 0 + 2*(m) + (n) )) 
#define Md_Om(m,n)   ( *(Hsys      + 4 + 2*(m) + (n) )) 
#define rho_in(m,n)  ( *(derivates + 0 + 2*(m) + (n) )) 
#define rho_out(m,n) ( *(derivates + 4 + 2*(m) + (n) )) 
 
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

  char date[256];  time_t curtime;   struct tm *loctime; curtime = time(NULL);  loctime = localtime (&curtime);
  char tme[1024];
  strftime (date, 256, "%d_%m", loctime); 
  strftime(tme,1024,"%H_%M_%S",loctime);
  
srand( time(NULL) );
double rn1 = (rand()%1000)*0.001;
double rn2 = (rand()%1000)*0.001;
double rn3 = (rand()%1000)*0.001;
double rn4 = (rand()%1000)*0.001;

printf("now: %d-%d-%d %d:%d:%d -- (%.3f,%.3f,%.3f,%.3f) \n", loctime->tm_year + 1900, loctime->tm_mon + 1, loctime->tm_mday, loctime->tm_hour, loctime->tm_min, loctime->tm_sec,rn1,rn2,rn3,rn4);
 
 
// ###########################################################################
// ########################## Coupling Element _ Rotating Frame LA ###########
// ###########################################################################

complex double *la_coup     = calloc(2*N_b,sizeof(double));
        double *n_ph        = calloc(  N_b,sizeof(double));
complex double *Hsys        = calloc(2*20,sizeof(double));
complex double *derivates   = calloc(2*20,sizeof(double));
		

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
printf("q_max=%.10f -- Polaronshift=%.10f \n",float_q,la_polaron);
#ifdef PHONON_COUPLING_OUTPUT
fclose(f_coup);
#endif

loctime = localtime (&curtime);
strftime(tme,1024,"%H_%M_%S",loctime);
printf("%s \n",tme);

FILE * f_dat;
snprintf(FILE_NAME,2048+1024,"%s_%s_%iK_BLOCH_PROPAGATOR_w_Nb=%i_OmL=%.3f_PHON=%.2f.dat",date,tme,(int)TEMPERATUR,N_b,OmL,PHON_FAK); 
f_dat=fopen(FILE_NAME,"w");

complex double IBM_arg = 0.;
int n,cut;	
int ehour = loctime -> tm_hour; int eminute = loctime -> tm_min; int esecond = loctime -> tm_sec;
int eyear = 1900+loctime -> tm_year; int emonth = loctime -> tm_mon; int eday  = loctime -> tm_mday;
double time_before,time_now; time_before =esecond + 60*eminute+60*60*ehour;

complex double phon_tab_n[N_b];
complex double phon_tab_m[N_b];
complex double phon_fak[N_b];
for (k=0;k<N_b;k++)
{
	float_q = q_0 + delta_q * k;
	phon_tab_n[k] = 1.;
	phon_tab_m[k] = 1.;
	phon_fak[k] = cexp(I*c_q*float_q*delta_t);
}

complex double tcexp;
complex double Snm=0.;    
double shift=0.;
complex double spec_fak = 0.;
int ANA_SWITCH;
int rr,rl; 
int ll,lr;

rho_out(0,0) = 1.; rho_out(1,0) = 0.;
rho_out(1,1) = 0.; rho_out(0,1) = 0.;

complex double num11,num01;
complex double num11_in,num01_in;
complex double puls;
double int_puls,tp;
int ts;
double OmR;

for(k=1; k<N; k++) // for increases the index after the loop is completed, k=0 in the first step
   {  
   t=delta_t*k;
   Snm = 0. + I*0.;
   // just for benchmarking the ibm
   IBM_arg     = 0. + I*0.;
   // calculate argument just one every time step k
   ANA_SWITCH = 1;
   shift = 0.;
   rr = round((rand()%100)*0.01);
   rl = round((rand()%100)*0.01);
   lr = round((rand()%100)*0.01);
   ll = round((rand()%100)*0.01);
   
   int_puls = 0.;
   for(ts=0;ts<20;ts++)
   {  
       tp = (k-1)*delta_t + (delta_t/20.)*ts;
	   //int_puls += (delta_t*OmL/20.) * (5.0*M_PI/WIDTH)*cexp(-1.0*( tp -t_NULL)*(tp-t_NULL)*1.386/(WIDTH*WIDTH)); 
	   int_puls  += (1.*OmL/20.)* (5.0*M_PI/WIDTH);
   }
   int_puls = OmL;
   OmR = creal(int_puls);
   
   M_Om(0,0) = cos(OmR*delta_t); 
   M_Om(1,1) = cos(OmR*delta_t); 
   M_Om(0,1) = -I*sin(OmR*delta_t); 
   M_Om(1,0) = -I*sin(OmR*delta_t); 
   
   Md_Om(0,0)=M_Om(0,0); Md_Om(1,1)=M_Om(1,1);
   Md_Om(1,0)=conj(M_Om(0,1)); Md_Om(0,1)=conj(M_Om(1,0));
   
/*   for(i=0; i<N_b; i++)
      {
       float_q  = q_0 + delta_q * i;
       spec_fak = delta_q*la_coup[i]*la_coup[i]/(2.*PI*PI*c_q*c_q); 	
	   // benchmark IBM argument
       IBM_arg +=  delta_q * (1./(2.*M_PI*M_PI)) * la_coup[i]  * conj(la_coup[i]) *   I*float_q*t/c_q; 
       IBM_arg += -delta_q * (1./(2.*M_PI*M_PI)) * la_coup[i]  * conj(la_coup[i]) * ( (2.*n_ph[i]+1.)*( 1.-cos(c_q*float_q*t) ) + I*sin(c_q*float_q*t) ) /(c_q*c_q);
	     shift +=  delta_q * (1./(2.*M_PI*M_PI)) * la_coup[i]  * conj(la_coup[i]) *   float_q/c_q;
       // now calculate it with path integral propagator for a time step k
	   // initialize n rotation 
	   phon_tab_n[i] = 1.;
	   if (k<180) cut=k;
	   else cut=180;
	   
       for(n=0;n<cut;n++)
        {
		// initialize m rotation 
        phon_tab_m[i] = 1.;
          for(m=0;m<n;m++)
	        {
            // phonon kernel time rotation
		    tcexp=phon_tab_n[i]*conj(phon_tab_m[i]);	
 		    if(n==m) 
		       {
               Snm += -spec_fak*(1.-creal(phon_fak[i]))*(2.*n_ph[i]+1.);
               Snm += -spec_fak*(-I*delta_t*c_q*float_q+I*cimag(phon_fak[i]));
		 	   }
		    else 
 		      {
   		      Snm += -2.*spec_fak*(1.-creal(phon_fak[i]))*( 1.*creal(tcexp)*(2.*n_ph[i]+1.));
			  Snm += -2.*spec_fak*(1.-creal(phon_fak[i]))*( -I*cimag(tcexp) );
              } 			
            // increase phonon_fak_m
			phon_tab_m[i] *= phon_fak[i];
		    } // loop m
        // after one m loop no more analytical ibm
		phon_tab_n[i] *= phon_fak[i];
		} // n loop
	  }// i loop	
   if (k>180) Snm += I*shift*(k-180)*delta_t; */
 
   rho_in(0,0) = rho_out(0,0);   rho_in(0,1) = rho_out(0,1);
   rho_in(0,1) = rho_out(0,1);   rho_in(1,1) = rho_out(1,1);
   // apply Rabi rotations from the left
   rho_out(0,0)= M_Om(0,0)*rho_in(0,0) + M_Om(0,1) * rho_in(1,0);
   rho_out(0,1)= M_Om(0,0)*rho_in(0,1) + M_Om(0,1) * rho_in(1,1);
   rho_out(1,0)= M_Om(1,0)*rho_in(0,0) + M_Om(1,1) * rho_in(1,0);
   rho_out(1,1)= M_Om(1,0)*rho_in(0,1) + M_Om(1,1) * rho_in(1,1);
 
   rho_in(0,0) = rho_out(0,0);   rho_in(0,1) = rho_out(0,1);
   rho_in(0,1) = rho_out(0,1);   rho_in(1,1) = rho_out(1,1);

   // apply Rabi rotations from the right    
   rho_out(0,0)= rho_in(0,0)*Md_Om(0,0) + rho_in(0,1) * Md_Om(1,0);
   rho_out(0,1)= rho_in(0,0)*Md_Om(0,1) + rho_in(0,1) * Md_Om(1,1);
   rho_out(1,0)= rho_in(1,0)*Md_Om(0,0) + rho_in(1,1) * Md_Om(1,0);
   rho_out(1,1)= rho_in(1,0)*Md_Om(0,1) + rho_in(1,1) * Md_Om(1,1);

   num11_in = num11; num01_in=num01;
   num11  = num11_in  -2.*delta_t*OmR* cimag(num01_in)/sqrt(2.);
   num01  = num01_in  + I*delta_t*OmR*(2.*   num11_in - 1.)/sqrt(2.);
 
//   fprintf(f_dat,"%.10f  \t %.10f \t %.10f \n",t/1000.,cabs(IBM_arg),cabs(Snm) );
   if (cabs(num11)>2.00) num11=0.;
   fprintf(f_dat,"%.10f  \t %.10f \t %.10f \t %.10f \n",t/1000.,creal(OmR),creal(rho_out(1,1)),creal(num11));
   
   if (k%SHOW_EVERY_STEP==0)
    {
     printf("STEP=%i ",k);
//     printf("RE_ANA=%.5f -- RE_EFF=%.5f -- DEV=%.5f --",creal(IBM_arg),creal(Snm),100.*(creal(IBM_arg)-creal(Snm))/creal(IBM_arg) );
//     printf("IM_ANA=%.5f -- IM_EFF=%.5f -- DEV=%.5f --",1000.*cimag(IBM_arg)/t,1000.*cimag(Snm)/t,100.*(cimag(IBM_arg)-cimag(Snm))/cimag(IBM_arg) );
       printf("rho_00=%.5f -- rho_11=%.5f -- puls=%.10f ",creal(rho_out(0,0)),creal(rho_out(1,1)),creal(puls));
//     printf("(%i,%i,%i,%i) -- ",rr,rl,lr,ll);
     curtime = time(NULL); loctime = localtime (&curtime); 
     ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
     time_now = esecond + 60*eminute+60*60*ehour;  
     printf("Time elapsed=%.2f sec.",time_now-time_before);
     time_before = time_now;
     printf("\n");
    } 
   } // k loop    
		 
fclose(f_dat);    
free(la_coup);
free(n_ph);
return 0;
}

