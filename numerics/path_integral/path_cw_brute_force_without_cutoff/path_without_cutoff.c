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

 double delta_t = 400.00; // in fs //((double)TIME_DYNAMICS_TO_NS) * 1000. * 1000. /((double)T_STEPS); 
 #define CUT_OFF_DENSITY ( 0.00000001 * 5.33 * 1.  )
 // to reproduce IBM set cutoff zero, otherwise coherence will vanish as with every time step the coherence
 // is spread over many states which are from a certain number on too small to be accounted for
 #define N_DELTA_T ( 1000 )
 #define MEMORY_LENGTH ( 10 )
 #define WITH_PATH_INTEGRAL_DYNAMICS
 // ############# ELECTRONIC PARAMETERS (start) ########################################################
 #define AW_00 (    1.0  )
 #define AW_01 (  I*0.5 * 0. )
 #define AW_10 ( -I*0.5 * 0. )
 #define AW_11 (    0.0  )
 // x 1/fs --> 1000 x 1/ps 
 #define GOmL      ( 0.00075184 * 0.0  ) 
 #define OmL      ( 0.00001000000 *1.0  ) 
 #define DELTA    ( 0.0135*0. )
 #define t_NULL   ( delta_t *10. )
 #define WIDTH    ( 1.05*1.*2500.0) 

 #define PHON_FAK (2.0 * 0.)
 // ############# PHONONIC PARAMETERS (start) ##########################################################
 //#define PHONON_COUPLING_OUTPUT
 #define TEMPERATUR (77.0001)//(2.01)//(1.7500100001) 
 #define BOSE(x) ( 1./( exp( (x) /(kB*TEMPERATUR)) - 1. ))
 // the number of modes is mostly determined by the length of the time trajectory - 100 for 50ps e.g.
 #define N_b ( 55000 ) //(2000)
 #define q_0  (0.0000000001)
 #define q_N (7.50000)
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
    #define rho_mat (33514170.0) // eVfs^2 /nm^5
    // kg = ... nutze Energie=Kraft x Weg , [kg] = 6.24151*10^3 ev fs^2/nm^5  
    // eV =1.6021766208(98)×10^19 J 
    // rho 5370 kg/m^3 --> 5370*6241,51 eV fs^2/nm^5 = 33514170 eV fs^2/nm^5
    #define D_c (-14.6*1.0) // eV
    #define D_v (-4.8*1.0)  // eV
 
// ############# PHONONIC PARAMETERS (end) ############################################################

// ############# PATH_INTEGRAL PARAMETERS (end) ###################################################################
#define M_Om(m,n,i)  ( *(Hsys + 4*(i) + 2*(m) + (n) ) ) 
#define path(n,p)     ( *(REL_STATES_INDEX   + (n)*(Nc+1) + (p) ) ) 
#define path_2(n,p)   ( *(REL_STATES_INDEX_2 + (n)*(Nc+1) + (p) ) ) 

/*
4^ 1: 4        -- 4^ 2: 16        -- 4^ 3:     64 -- 4^ 4:     256 -- 4^ 5:    1024 -- 4^ 6:     4096 -- 
4^ 7: 16384    -- 4^ 8: 65536     -- 4^ 9: 262144 -- 4^10: 1048576 -- 4^11: 4194304 -- 4^12: 16777216 
4^13: 67108864 -- 4^14: 268435456 -- 4^15: 1073741824 */

		
// #######################################################
// ############# PHONON COUPLING ELEMENTS ################
// #######################################################

double gvc (double q)
 { double faktor,argument_h,argument_el;
   double temp_gvc; 
   argument_h  = -1.* hbar * q * hbar * q * 0.25 / (hbar_omega_h * eff_mass_h   * el_mass );
   argument_el = -1.* hbar * q * hbar * q * 0.25 / (hbar_omega_el * eff_mass_el * el_mass );
   faktor   = sqrt (0.5 * hbar * q / ( rho_mat * c_q) );
   temp_gvc =  faktor * ( D_c * exp( argument_el) - D_v * exp( argument_h) ) * PHON_FAK;
   return temp_gvc;  
 } 


void INITIALIZE_PHONON_INTERACTION(complex double *la_coup,double *n_ph)
{
double la_polaron = 0.;
double float_q = 0.;
int i;

for (i=0;i<N_b;i++)
{
  float_q      = q_0 + delta_q * i;
  la_coup[i]   = gvc(float_q);
  la_polaron  += (delta_q/(2.*PI*PI))*(float_q/c_q) * la_coup[i] * conj(la_coup[i]);
  n_ph[i]      = BOSE((hbar*c_q*float_q)); 
}
printf("q_max=%.10f -- Polaronshift=%.10f \n",float_q,la_polaron);

#ifdef PHONON_COUPLING_OUTPUT
FILE * f_coup;
snprintf(FILE_NAME,2048+1024,"%s_%iK_PAR_LA_phon_COUPLING_w_Nb=%i_OmL=%.3f_PHON=%.2f.dat",date,(int)TEMPERATUR,N_b,OmL,PHON_FAK); 
f_coup=fopen(FILE_NAME,"w");
for (i=0;i<N_b;i++) fprintf(f_coup,"%.10f \t %.10f \t %.10f \n",q_0 + delta_q * i,cabs(la_coup[i]),n_ph[i]);
fclose(f_coup);
#endif

return;
}

// #######################################################
// #######################################################
// #######################################################

void PRINT_PATH(int *PATH, int Nc)
{
int i; 
printf("(%*i",2,PATH[0]);  for (i=1;i<Nc+1;i++) printf(",%*i",2,PATH[i]); printf(")\n"); 
printf("(%*i",2,0);  for (i=1;i<Nc+1;i++) printf(",%*i",2,i); printf(")\n"); 
}

// #######################################################
// #######################################################
// #######################################################


// #######################################################
// ################## MAIN ###############################
// #######################################################

int main(){  
int k,i,m;
double t;

char date[256];  time_t curtime;   struct tm *loctime; curtime = time(NULL);  loctime = localtime (&curtime);
char tme[1024]; strftime (date, 256, "%d_%m", loctime); strftime(tme,1024,"%H_%M_%S",loctime);
  
srand( time(NULL) );
double rn1 = (rand()%1000)*0.001;
int hour=loctime->tm_hour;
int sec =loctime->tm_sec ;
int min =loctime->tm_min ;
printf("now: %d-%d-%d %d:%d:%d -- (Random %.3f) \n", loctime->tm_year + 1900, loctime->tm_mon + 1, loctime->tm_mday, loctime->tm_hour, loctime->tm_min, loctime->tm_sec,rn1);
char FILE_NAME[2048+1024];
 
// ###########################################################################
// ########################## PATH_MAPPING FOR AUGMENTED DENSITY MATRIX ######
// ###########################################################################
// there are four options per time increment 0: 0 - l0,r0;  1 - l0,r1;  2 - l1,r0;  3 - l1,r1; 

int Nc=MEMORY_LENGTH; 
int PATH[Nc+2]; 

// ###########################################################
// ############ Define Bloch rotations #######################
        double *PULSE_AREA  = calloc(  (N_DELTA_T+2*Nc)  , sizeof(double));
        double *PULSE       = calloc(  (N_DELTA_T+2*Nc)  , sizeof(double));
complex double *Hsys        = calloc(2*(N_DELTA_T+2*Nc)*4, sizeof(double));

double OmR;
double int_puls,tp,puls_max;
puls_max = -1.;
int ts;
for(k=0;k<(N_DELTA_T+2*Nc);k++)
{	
int_puls = 0.;
for(ts=0;ts<20;ts++)
   {  
   tp = (k-1)*delta_t + (delta_t/Nc)*ts;
   int_puls += (delta_t/Nc) * (5.0*GOmL/WIDTH)*cexp(-1.0*( tp-t_NULL)*(tp-t_NULL)*1.386/(WIDTH*WIDTH)); 
   }
//printf("%.10f -- %.10f   \n",k*delta_t*0.001,int_puls); 
OmR = creal(int_puls);
OmR = creal(OmL)     ; // ATTENTION now cw pulse // comment to go to Gaussian
PULSE[k] = creal(int_puls) ;
if (PULSE[k]>puls_max) puls_max = PULSE[k];
if (k==0) PULSE_AREA[k] = delta_t * OmR;
else PULSE_AREA[k] = delta_t * OmR + PULSE_AREA[k-1];
//OmR = delta_t*OmL;
M_Om(0,0,k) =    cos(OmR*delta_t);
M_Om(1,1,k) =    cos(OmR*delta_t);
M_Om(0,1,k) = -I*sin(OmR*delta_t);
M_Om(1,0,k) = -I*sin(OmR*delta_t);  
}   
// ############ End Define Bloch rotations #######################
// ###############################################################

// ###########################################################
// ############ Define Phonon Kernel #########################
complex double *la_coup     = calloc(2*N_b,sizeof(double));
        double *n_ph        = calloc(  N_b,sizeof(double));
complex double PHONON_KERNEL[Nc+10];
int p;
double float_q,spec_fak,IBM_arg;

INITIALIZE_PHONON_INTERACTION(la_coup,n_ph);

// ---------------------------------------
// first calculate as a benchmark the IBM
complex double IBM[N_DELTA_T+2*Nc];
IBM[0] = cimag(AW_01);
for(p=1;p<(N_DELTA_T+2*Nc);p++)
{
IBM_arg = 0.;
for(i=0; i<N_b; i++)
      {
       float_q  = q_0 + delta_q * i;
       IBM_arg += -delta_q * (1./(2.*M_PI*M_PI)) * la_coup[i] * conj(la_coup[i]) * (2.*n_ph[i]+1.)*( 1.-cos(c_q*float_q*p*delta_t) ) /(c_q*c_q);  
      }
IBM[p] = cimag(AW_01)*cexp(IBM_arg);
}
// end IBM -------------------------------
// ---------------------------------------

for(p=0;p<Nc+10;p++)
{
PHONON_KERNEL[p] = 0. +I*0.;	
for(i=0; i<N_b; i++)
      {
       float_q  = q_0 + delta_q * i;
       spec_fak = delta_q * (2./(2.*M_PI*M_PI*c_q*c_q)) * la_coup[i]  * conj(la_coup[i]) ; 	
	   if (p==0)
   	   PHONON_KERNEL[p] += -    0.5*spec_fak* (2.*n_ph[i]+1.) * ( 1.-cos(c_q*float_q*delta_t) ) 
                           - I* 0.5*spec_fak*                        sin(c_q*float_q*delta_t)  ;
	   else 
	   PHONON_KERNEL[p] += -  1.0*spec_fak* (2.*n_ph[i]+1.) * ( 1.-cos(c_q*float_q*delta_t) ) * cos(c_q*float_q*delta_t*p) 
                           +I*1.0*spec_fak*                   ( 1.-cos(c_q*float_q*delta_t) ) * sin(c_q*float_q*delta_t*p)  ;    
      }
printf("PHONON_KERNEL[%i]=%.15f -- IBM=%.10f \n",p,creal(PHONON_KERNEL[p]),creal(IBM[p]) );
}

printf("dt=%.10ffs - dt=%.10fps -- %.10f \n",delta_t,delta_t*0.001,fabs(PHONON_KERNEL[20]) );
if (cabs(PHONON_KERNEL[20])>0.000001) printf("ATTENTION!!! MEMORY LENGTH MIGHT BE INSUFFICIENT!!! \n");

// ############ End Phonon Kernel ############################
// ###########################################################

#ifdef WITH_PATH_INTEGRAL_DYNAMICS

double *DYNAMICS_01 = calloc( (N_DELTA_T+2*Nc) , sizeof(double));
double *DYNAMICS_11 = calloc( (N_DELTA_T+2*Nc) , sizeof(double));
double *RE_POL      = calloc( (N_DELTA_T+2*Nc) , sizeof(double));

int MAX_NUMBER_OF_STATES = 419430;

int *REL_STATES_INDEX   = calloc( (Nc+1)*MAX_NUMBER_OF_STATES, sizeof(int));
int *REL_STATES_INDEX_2 = calloc( (Nc+1)*MAX_NUMBER_OF_STATES, sizeof(int));

//complex double init_rho;
complex double test=1.;
//int NEGLIGIBLE =0;
int NON_NEGLIGIBLE = 0;
//double OmRp;
//int ind0,ind1,L_LIMIT,K_LIMIT;
int MAX_NN,l;
int Dt = N_DELTA_T;


// ###########################################################################
// ###########################################################################
// ###########################################################################
// ###########################################################################
// ###########################################################################

complex double PATH_INIT[4]; // is the initial state for given calculation
complex double PATH_INIT_2[4]; // saves the next initial state
complex double RESULT[4]; // captures the given density matrix result
int HELPER = 0; // to switch between relevant path arrays

PATH_INIT_2[0]=AW_00;   
PATH_INIT_2[1]=AW_01;   
PATH_INIT_2[2]=AW_10;   
PATH_INIT_2[3]=AW_11;   

DYNAMICS_01[0]=cimag(PATH_INIT[1]);
DYNAMICS_11[0]=creal(PATH_INIT[2]);

for(p=0;p<4;p++) printf("[%i] --> (%i,%i) \n",p,p/2,p%2);

printf("---------------------------------------------------------------------- \n");
printf("------- NOW, EVERY TIME STEP INCLUDES FULL MEMORY LOOP  -------------- \n");
printf("---------------------------------------------------------------------- \n");

// ####################################################
// ####################################################
// ########### TIME INTEGRATION #######################
// ####################################################
// ####################################################
double kernel_temp = 0.;
for(i=0; i<N_b; i++)
      {
       float_q  = q_0 + delta_q * i;
       spec_fak = -1. * delta_q * (1./(2.*M_PI*M_PI*c_q*c_q)) * la_coup[i]  * conj(la_coup[i]) * (2.*n_ph[i]+1.) ;  
       kernel_temp += spec_fak * ( cos(c_q*float_q*  Nc     *delta_t) - cos(c_q*float_q*(Nc-1.)*delta_t)  ) ;  
       kernel_temp += spec_fak * ( cos(c_q*float_q*( 1.-1. )*delta_t) - cos(c_q*float_q* 1.    *delta_t)  ) ;          
      }

complex double Snm = 1.;
int DEPTH = 0;

for(Dt=0;Dt<N_DELTA_T;Dt++)
{
if (Dt>0) PHONON_KERNEL[19] -= kernel_temp;    
// the same again with new initial state
for(l=0;l<4;l++) PATH_INIT[l]=PATH_INIT_2[l];
//printf("00=%.10f \n",creal(PATH_INIT[0]));printf("01=%.10f \n",cimag(PATH_INIT[1]));printf("10=%.10f \n",cimag(PATH_INIT[2]));printf("11=%.10f \n",creal(PATH_INIT[3]));
// to include four initial states
NON_NEGLIGIBLE = 4; 
path(0,0)=0;   path(1,0)=1;   path(2,0)=2;   path(3,0)=3;
HELPER=0;
// always starting at p=1, because PATH[0] is the initial state!!!!
for(p=1;p<Nc;p++)
{
MAX_NN = NON_NEGLIGIBLE;
NON_NEGLIGIBLE = 0;
for(l=0;l<4;l++) RESULT[l]=0.;

for (k=0;k<MAX_NN;k++)
{	
 // first write path array   
 if (HELPER==0) { for(m=0;m<Nc;m++) PATH[m]=  path(k,m); }
 else           { for(m=0;m<Nc;m++) PATH[m]=path_2(k,m); }
 // now calculate the next time step 
 for(l=0;l<4;l++)
 {
 PATH[p] = l ;   
 test = 1.; Snm = 1.;
 for(i=0;i<p;i++)
 {    
  // include the initial state of the adm
    if (i==0)  test *= PATH_INIT[PATH[0]];   
  // calculate the next system time step      
    test *= M_Om((PATH[i+1]/2),(PATH[i]/2),i) * conj(M_Om((PATH[i+1]%2),(PATH[i]%2),i));
    for(m=0;m<=i;m++) { Snm *= cexp( ((PATH[i+1]/2) - (PATH[i+1]%2)) * ((PATH[i+1-m]/2)*PHONON_KERNEL[m]-(PATH[i+1-m]%2)*conj(PHONON_KERNEL[m])) ); } 
 }
   // for(m=0;m<=p;m++) { Snm *= cexp( ((PATH[p]/2) - (PATH[i+1]%2)) * ((PATH[i+1-m]/2)*PHONON_KERNEL[m]-(PATH[i+1-m]%2)*conj(PHONON_KERNEL[m])) ); } 
 
 if ( ( cabs(test) > CUT_OFF_DENSITY ) )   
    {
         RESULT[l] += test * Snm; //if (l==1) DEPTH++;
         // store the path 
         //printf("PATH[%i]:( ",NON_NEGLIGIBLE-1); for(m=0;m<Nc;m++) printf("%i,",path(k,m)); printf(")\n");
		 if (HELPER == 0) { for(m=0;m<Nc;m++) path_2(NON_NEGLIGIBLE,m)=PATH[m]; }
		 else             { for(m=0;m<Nc;m++)   path(NON_NEGLIGIBLE,m)=PATH[m]; }
         //printf("PATH[%i]:( ",NON_NEGLIGIBLE); for(m=0;m<Nc;m++) printf("%i,",path_2(NON_NEGLIGIBLE,m)); printf(")\n");
         NON_NEGLIGIBLE++;
     }
} // l loop four states 0,1,2,3 for new time step
} // k loop all relevant states

if (HELPER==0) HELPER=1;
else HELPER = 0;

if (p==1) {for(l=0;l<4;l++) PATH_INIT_2[l]=RESULT[l];}
if (Dt==0)
{
m=p;
printf("[%i]:t=%.5fps: ADM_00=%.10f -- BLOCH=%.10f -- DIFF=%.10f\n",m,m*delta_t*0.001,creal(RESULT[0]),cos(PULSE_AREA[m-1])*cos(PULSE_AREA[m-1]),creal(RESULT[0])-cos(PULSE_AREA[m-1])*cos(PULSE_AREA[m-1]) );
printf("[%i]:t=%.5fps: RE[01]=%.10f -- BLOCH=%.10f -- DIFF=%.10f \n",m,m*delta_t*0.001,creal(RESULT[1]),0.5*sin(2.*PULSE_AREA[m-1]),cimag(RESULT[1])-0.5*sin(2.*PULSE_AREA[m-1]) );
printf("[%i]:t=%.5fps: IM[01]=%.10f -- IBM  =%.10f -- DIFF=%.10f \n",m,m*delta_t*0.001,cimag(RESULT[1]),creal(IBM[p-1]),creal(IBM[p-1])-cimag(RESULT[1]) );
printf("[%i]:t=%.5fps: ADM_11=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",m,m*delta_t*0.001,creal(RESULT[3]),sin(PULSE_AREA[m-1])*sin(PULSE_AREA[m-1]),creal(RESULT[3])-sin(PULSE_AREA[m-1])*sin(PULSE_AREA[m-1]) );
printf("[%i]:t=%.5fps: PULSE =%.10f -- DEPTH=%i \n\n",m,m*delta_t*0.001,creal(PULSE[m-1]),DEPTH);
DYNAMICS_01[p]=cimag(RESULT[1]);
DYNAMICS_11[p]=creal(RESULT[3]);
     RE_POL[p]=creal(RESULT[1]);
printf("RELEVANT STATES %i \n",NON_NEGLIGIBLE);
}
} // p loop

if(Dt>0)
{    
m=Dt+p-1;
printf("[%i]:t=%.5fps: ADM_00=%.10f -- BLOCH=%.10f -- DIFF=%.10f\n",m,m*delta_t*0.001,creal(RESULT[0]),cos(PULSE_AREA[m-1])*cos(PULSE_AREA[m-1]),creal(RESULT[0])-cos(PULSE_AREA[m-1])*cos(PULSE_AREA[m-1]) );
printf("[%i]:t=%.5fps: RE[01]=%.10f -- BLOCH=%.10f -- DIFF=%.10f \n",m,m*delta_t*0.001,creal(RESULT[1]),0.5*sin(2.*PULSE_AREA[m-1]),cimag(RESULT[1])-0.5*sin(2.*PULSE_AREA[m-1]) );
printf("[%i]:t=%.5fps: IM[01]=%.10f -- IBM  =%.10f -- DIFF=%.10f \n",m,m*delta_t*0.001,cimag(RESULT[1]),creal(IBM[p-1]),creal(IBM[p-1])-cimag(RESULT[1]) );
printf("[%i]:t=%.5fps: ADM_11=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",m,m*delta_t*0.001,creal(RESULT[3]),sin(PULSE_AREA[m-1])*sin(PULSE_AREA[m-1]),creal(RESULT[3])-sin(PULSE_AREA[m-1])*sin(PULSE_AREA[m-1]) );
printf("[%i]:t=%.5fps: PULSE =%.10f -- DEPTH=%i \n\n",m,m*delta_t*0.001,creal(PULSE[m-1]),DEPTH);
DYNAMICS_01[m]=cimag(RESULT[1]);
DYNAMICS_11[m]=creal(RESULT[3]);
     RE_POL[m]=creal(RESULT[1]);
printf("RELEVANT STATES %i \n",NON_NEGLIGIBLE);
}
} // Dt loop 

// ####################################################
// ####################################################
// ########### END OF TIME INTEGRATION ################
// ####################################################
// ####################################################



free(REL_STATES_INDEX);
free(REL_STATES_INDEX_2);


snprintf(FILE_NAME,2048+1024,"%d_Bloch_Path_Integral_t=%.1fps_CutOff=%.5f_Nc=%i_DETUNING=%.6f_OmL=%.5f_GOmL=%.5f_%.0fK.dat",100*100*hour+100*min+sec,m*delta_t*0.001,CUT_OFF_DENSITY/(GOmL+OmL),Nc,DELTA/(GOmL+OmL),OmL,GOmL,TEMPERATUR); 
FILE *f_dat; f_dat = fopen(FILE_NAME,"w");
for(k=0;k<Nc+Dt;k++)
{
fprintf(f_dat,"%.5f \t",delta_t*k*0.001);
fprintf(f_dat,"%.20f  \t %.20f  \t %.20f \t",RE_POL[k],DYNAMICS_01[k],DYNAMICS_11[k]);
fprintf(f_dat,"%.10f \t %.10f \t",sin(PULSE_AREA[k-1])*sin(PULSE_AREA[k-1]),0.5*sin(2.*PULSE_AREA[k-1]));
fprintf(f_dat,"%.10f \t %.10f \t",sin(PULSE_AREA[k-1])*sin(PULSE_AREA[k-1])-DYNAMICS_11[k],0.5*sin(2.*PULSE_AREA[k-1])-DYNAMICS_01[k]);
//fprintf(f_dat,"%.5f",PULSE[k]/puls_max);
fprintf(f_dat,"\n");
}
fclose(f_dat); 

printf("dt=%.10ffs - dt=%.10fps \n",delta_t,delta_t*0.001 );

free(DYNAMICS_01); 
free(DYNAMICS_11); 
free(RE_POL);

#endif

free(PULSE_AREA);
free(la_coup);
free(n_ph);
return 0;
}

