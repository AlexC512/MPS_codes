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
 #define T_STEPS 50
 int N=T_STEPS;  
 double delta_t = ((double)TIME_DYNAMICS_TO_NS) * 1000. * 1000. /((double)T_STEPS); 
 #define MEMORY_CUT (20.0*1000.0)
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
 //#define PHONON_COUPLING_OUTPUT
 #define TEMPERATUR (4.0001)//(2.01)//(1.7500100001) 
 #define BOSE(x) ( 1./( exp( (x) /(kB*TEMPERATUR)) - 1. ))
 // the number of modes is mostly determined by the length of the time trajectory - 100 for 50ps e.g.
 #define N_b ( 100 ) //(2000)
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

// #######################################################
// #######################################################
// ####################################################### 
 complex double IBM_BRUTE_FORCE
 (int k, 
  int cut,
  int i, 
  complex double phon_fak, 
          double n_ph,   
  complex double spec_fak)
 {
 complex double Snm_i = 0.+I*0.;
 int m,n; 
 double float_q = q_0 + delta_q * i;
 
 complex double phon_tab_n = cexp(I*c_q*float_q*delta_t*cut);
 complex double phon_tab_m = 1.+0.*I;
 complex double phon_rot;

       for(n=cut;n<k+1;n++)
	   {
		phon_tab_m = 1.+0.*I;
          for(m=0;m<n+1;m++)
	        {
            // phonon kernel time rotation
		    phon_rot = phon_tab_n*conj(phon_tab_m);	
 		    if(n==m) 
		       {
               Snm_i += -spec_fak*(1.-creal(phon_fak))*(2.*n_ph+1.);
               Snm_i += -spec_fak*(-I*delta_t*c_q*float_q+I*cimag(phon_fak));
		 	   }
		    else 
 		      {
   		      Snm_i += -2.*spec_fak*(1.-creal(phon_fak))*( 1.*creal(phon_rot)*(2.*n_ph+1.));
			  Snm_i += -2.*spec_fak*(1.-creal(phon_fak))*( -I*cimag(phon_rot) );
              } 			
            // increase phonon_fak_m
			phon_tab_m *= phon_fak;
		    } // loop m
        // after one m loop no more analytical ibm
		phon_tab_n *= phon_fak;
		} // n loop */
 
 return Snm_i;	 
 }	 
// #######################################################
// #######################################################
// #######################################################

// #######################################################
// #######################################################
// #######################################################
 complex double IBM_EFFICIENT(int k, int i, complex double phon_fak, double n_ph, complex double spec_fak)
 {
 complex double Snm_i = 0.+I*0.;
 int m; 
 double float_q = q_0 + delta_q * i;
 
 complex double phon_tab_k = cexp(I*delta_t*k*c_q*float_q);
 complex double phon_tab_m = 1.+0.*I;
 complex double phon_rot;

 for(m=0;m<k+1;m++)
    {
     // phonon kernel time rotation
	 phon_rot = phon_tab_k*conj(phon_tab_m);	
 	 if(k==m) 
	 {
      Snm_i += -spec_fak*(1.-creal(phon_fak))*(2.*n_ph+1.);
      Snm_i += -spec_fak*(-I*delta_t*c_q*float_q+I*cimag(phon_fak));
     }
	 else 
 	 {
   	  Snm_i += -2.*spec_fak*(1.-creal(phon_fak))*( 1.*creal(phon_rot)*(2.*n_ph+1.));
	  Snm_i += -2.*spec_fak*(1.-creal(phon_fak))*( -I*cimag(phon_rot) );
     } 			
     // increase phonon_fak_m
	phon_tab_m *= phon_fak;
    } // loop m
 
 return Snm_i;	 
 }	 
// #######################################################
// #######################################################
// #######################################################


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
printf("now: %d-%d-%d %d:%d:%d -- (Random %.3f) \n", loctime->tm_year + 1900, loctime->tm_mon + 1, loctime->tm_mday, loctime->tm_hour, loctime->tm_min, loctime->tm_sec,rn1);
char FILE_NAME[2048+1024];
 
// ###########################################################################
// ########################## Coupling Element _ Rotating Frame LA ###########
// ###########################################################################
complex double *la_coup     = calloc(2*N_b,sizeof(double));
        double *n_ph        = calloc(  N_b,sizeof(double));

double la_polaron = 0.;
double float_q = 0.;

for (i=0;i<N_b;i++)
{
  float_q      =q_0 + delta_q * i;
  la_coup[i]   = gvc(float_q);
  la_polaron  += (delta_q/(2.*PI*PI))*(float_q/c_q) * la_coup[i] * conj(la_coup[i]);
  n_ph[i] = BOSE((hbar*c_q*float_q)); 
}
printf("q_max=%.10f -- Polaronshift=%.10f \n",float_q,la_polaron);

#ifdef PHONON_COUPLING_OUTPUT
FILE * f_coup;
snprintf(FILE_NAME,2048+1024,"%s_%iK_PAR_LA_phon_COUPLING_w_Nb=%i_OmL=%.3f_PHON=%.2f.dat",date,(int)TEMPERATUR,N_b,OmL,PHON_FAK); 
f_coup=fopen(FILE_NAME,"w");
for (i=0;i<N_b;i++) fprintf(f_coup,"%.10f \t %.10f \t %.10f \n",q_0 + delta_q * i,cabs(la_coup[i]),n_ph[i]);
fclose(f_coup);
#endif

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

// ###########################################################################
// ########################## Coupling Element _ Rotating Frame LA ###########
// ###########################################################################

// ###########################################################################
// ########################## Dynamics #######################################
// ###########################################################################

loctime = localtime (&curtime);
strftime(tme,1024,"%H_%M_%S",loctime);
printf("Now: %s \n",tme);

FILE * f_dat;
snprintf(FILE_NAME,2048+1024,"%s_%s_%iK_IBM_Benchmark_PROPAGATOR_w_Nb=%i_OmL=%.3f_PHON=%.2f.dat",date,tme,(int)TEMPERATUR,N_b,OmL,PHON_FAK); 
f_dat=fopen(FILE_NAME,"w");


complex double tcexp;
complex double Snm=0.;    
double shift=0.;
complex double spec_fak = 0.;
int ANA_SWITCH;

complex double IBM_prev=0.;
double mem_depth = -10.9;
double steady_state;

int CUT_OFF = round( MEMORY_CUT/delta_t ); 
printf("Nc = %i \n",CUT_OFF);

complex double num_pol = 1.;
complex double ana_pol = 0.;

complex double IBM_arg = 0.;
int n,cut;	
int ehour = loctime -> tm_hour; int eminute = loctime -> tm_min; int esecond = loctime -> tm_sec;
int eyear = 1900+loctime -> tm_year; int emonth = loctime -> tm_mon; int eday  = loctime -> tm_mday;
double time_before,time_now; time_before =esecond + 60*eminute+60*60*ehour;

Snm = 0. + I*0.;
 
 
// there are four options per time increment 0: 0 - l0,r0;  1 - l0,r1;  2 - l1,r0;  3 - l1,r1; 
// [0,4,16] --> all ln=0,rn=0
// [1] --> l1=0,r1=1, all ln,rn=0; [3] --> l1=1,r1=1, all ln,rn=0 n>1
// [9] --> l1=0,r1=1,

// mapping with 4 powers,
int Nc=14; 
int PATH1[Nc]; 
PATH1[0] = 0; // we fix the initial conditions rho(t=0)= |0><0|, therefore paths start k=4
long int kl,k_map;
long int npath=0;
long int pow4[Nc+1]; pow4[0]=1;
for(i=1;i<Nc+1;i++) { pow4[i] = 4*pow4[i-1]; printf("4^%*i: %i \n",2,i,pow4[i]); }  
// the algorithm which maps an integer to a state path (0,1,3,0,2,3,1,...)
// the problem is, that even a long double allows only with 4^n a length of 14 (4^15-1) which is too short
// so I have 4^(Nc) - 1 paths consistent with the initial condition, maximum of Nc = 14, as I need 
// 4 powers up to 4^15 to map them. Going beyond, will not help, I need long int. 

// store the k numbers of the p, which are consistent with initial conditions
int *KPC  = calloc( (pow4[Nc-1]-1),sizeof(long int));

 for (kl=1;kl<pow4[Nc];kl++)
{	// for every k I need a corresponding path
	k_map = kl; 
	PATH1[0]=0;
	// check first, whether number is consistent with CHOSEN initial condition
	if ( k_map%4 == PATH1[0])
	{
	 KPC[npath]=kl;
     npath++; 	
	 k_map -= PATH1[0];
	 i = Nc;
	 while( (k_map>0) && (i>0) )  
	 { 
        if ( pow4[i] > k_map ) PATH1[i]=0; 
		else 			
		{
	     PATH1[i] = k_map/pow4[i] ;
		 k_map -= pow4[i]*PATH1[i];
     	}
	 i--; 	
     }
	 while ( i>0 ) { PATH1[i]=0; i--;}
	//check the mapping
	/*if (kl%4096==0)	{	
	printf("%*i ---------------------------------\n",3,npath);
	printf("%*i:(%i",4,kl,PATH1[0]);
	for (m=1;m<Nc;m++) printf(",%i",PATH1[m]);
	printf(")\n");
	}
	// */
    }
	//else printf("k=%i inconsistent with initial conditions \n",k);
}
printf("%*i:(%i",4,kl,PATH1[0]);
for (m=1;m<Nc;m++) printf(",%i",PATH1[m]);
printf(")\n");

printf("npath=%i -- \n",npath);
printf("     =%i -- \n",pow4[Nc-1]-1);

printf("OK, let's check the paths \n");
for (kl=1;kl<pow4[Nc-1]-1;kl++)
{	k_map = KPC[kl];
	if ( k_map%4 > 0 ) printf("k=%i inconsistent with initial conditions \n",k);
}
printf("Checking complete.\n");
// meaning, I have now all the integeres consistent with the initial conditions.
// to map path length greater than 14 which unfortunately might be necessary, I need an array
// for i to max_path, for k to max_path, map(i,k) which gives me a length of 28.
// I need this array to store my density matrices
// every rho(i,k) has a well-defined and associated state number( or path) assigned to it.


/* 
for(k=0; k<N; k++) // for increases the index after the loop is completed, k=0 in the first step
   {  
   t=delta_t*k;
   // just for benchmarking the ibm
   IBM_prev    = IBM_arg ;
   IBM_arg     = 0. + I*0.;
//   Snm = 0. + I*0.;
   // calculate argument just one every time step k
   shift = 0.;
   for(i=0; i<N_b; i++)
      {
       float_q  = q_0 + delta_q * i;
       spec_fak = delta_q*la_coup[i]*la_coup[i]/(2.*PI*PI*c_q*c_q); 	
	   // benchmark IBM argument
       IBM_arg +=  delta_q * (1./(2.*M_PI*M_PI)) * la_coup[i]  * conj(la_coup[i]) *   I*float_q*t/c_q; 
       IBM_arg += -delta_q * (1./(2.*M_PI*M_PI)) * la_coup[i]  * conj(la_coup[i]) * ( (2.*n_ph[i]+1.)*( 1.-cos(c_q*float_q*t) ) + I*sin(c_q*float_q*t) ) /(c_q*c_q);
	     shift +=  delta_q * (1./(2.*M_PI*M_PI)) * la_coup[i]  * conj(la_coup[i]) *   float_q/c_q;
       // now calculate it with path integral propagator for a time step k
	   if (k>CUT_OFF) cut=k-CUT_OFF;
	   else cut=0;
   	   //Snm += IBM_EFFICIENT(k,i,phon_fak[i],n_ph[i],spec_fak);
	   if (i==0) Snm=0.; Snm += IBM_BRUTE_FORCE(k,cut,i,phon_fak[i],n_ph[i],spec_fak);
	  }// i loop	
   if (cut>0) { Snm = I*shift*k*delta_t;  }
 
   if (k>0) fprintf(f_dat,"%.10f  \t %.10f \t %.10f \n",t/1000.,creal(IBM_arg),creal(Snm) );

   steady_state = 100.*(creal(IBM_arg) - creal(IBM_prev))/creal(IBM_arg);   
   if ( ( cabs(steady_state)<0.01000) && (k>50) && (mem_depth<0.1) )  mem_depth=k*delta_t;  

   num_pol  = cexp(Snm);
   ana_pol  = cexp(IBM_arg);

   if (k%1==0)
    {
     printf("%i:t=%.2f - ",k,k*delta_t);
     printf("RE_ANA=%.7f -- RE_EFF=%.7f -- DEV=%.5f --",creal(IBM_arg),creal(Snm),100.*(creal(IBM_arg)-creal(Snm))/creal(IBM_arg) );
     printf("IM_ANA=%.7f -- IM_EFF=%.7f -- DEV=%.5f --",cimag(IBM_arg),cimag(Snm),100.*(cimag(IBM_arg)-cimag(Snm))/cimag(IBM_arg) );
//     printf("POL_AN=%.7f -- POL_NM=%.7f -- MEM_DEPTH=%.5fps --",cabs(ana_pol),cabs(num_pol),mem_depth/1000.);
	 
     //printf("(%i,%i,%i,%i) -- ",rr,rl,lr,ll);
     curtime = time(NULL); loctime = localtime (&curtime); 
     ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
     time_now = esecond + 60*eminute+60*60*ehour;  
     printf("(%.2f).",time_now-time_before);
     time_before = time_now;
     printf("\n");
    }
   } // k loop    
  */	
		 
fclose(f_dat);    
free(la_coup);
free(KPC);
free(n_ph);
return 0;
}

