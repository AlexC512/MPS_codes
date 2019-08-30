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
 #define T_STEPS 250
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
 
 #define OmL      ( 0.0000784 * 2.0  ) 
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

#define M_Om(m,n)    ( *(Hsys      + 0 + 2*(m) + (n) )) 
#define Md_Om(m,n)   ( *(Hsys      + 4 + 2*(m) + (n) )) 

#define PATHS(m,n) ( *(path_table + (m)*(Nc+1) + (n) ) )
		 
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

void INITIALIZE_PHONON_INTERACTION(complex double *la_coup,double *n_ph, complex double *phon_fak)
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
  phon_fak[i]  = cexp(I*c_q*float_q*delta_t);
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

// #######################################################
// #######################################################
// #######################################################
void CREATE_PATH_TABLE(int *path_table,int MAX_K,int Nc,int *pow4,int INIT_RHO)
{

int PATH1[Nc]; 
PATH1[0] = INIT_RHO; // we fix the initial conditions rho(t=0)= |0><0|, therefore paths start k=4
int i,m,kl,k_map;
int npath=0;
//PATHS(npath,1)=PATH1[0];
//for(m=1;m<Nc;m++) PATHS(npath,m+1)=0;
//npath++;
int Nz=1; // I have to include before every new power of 4, an empty row
// to take into account that the groundstate is created several times
for (kl=0;kl<MAX_K;kl++)
{	// for every k I need a corresponding path
	k_map = kl; 
	PATH1[0]=INIT_RHO; 
	// check first, whether number is consistent with CHOSEN initial condition
	if ( k_map%4 == PATH1[0])
	{
	 PATHS(npath,0) = kl;  	
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
	 // here I fill in the zeros 
	 while ( i>0 ) { PATH1[i]=0; i--;}
	 // here before any powermanifold an empty row because of ground state calculations
	 // remember to use a path not consistent with initial conditions!!! but I stay with rho_00(0) anyhow
	 if ( (kl%pow4[Nz]==0) && (kl>0) ) { Nz++; for(m=1;m<Nc+1;m++) PATHS(npath,m)=0.; npath++; } 
     PATHS(npath,0)=kl+1;
	 for(m=0;m<Nc;m++) PATHS(npath,m+1)=PATH1[m];
	 npath++;
	}
}
/*
printf("END of MAPPING\n");
for(kl=0;kl<npath;kl++) 
{
 //if (kl%1==0)    { 
 	 printf("%*i:(%i",4,PATHS(kl,0),PATHS(kl,1)); 
	 for (m=2;m<Nc+1;m++) printf(",%i",PATHS(kl,m)); 
	 printf(")\n"); 
    //}	
}
printf("Number of possible path in one array\nnpath=%i -- \n",npath);
printf("     =%i -- \n",pow4[Nc-1]-1);
//*/

return;
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
printf("now: %d-%d-%d %d:%d:%d -- (Random %.3f) \n", loctime->tm_year + 1900, loctime->tm_mon + 1, loctime->tm_mday, loctime->tm_hour, loctime->tm_min, loctime->tm_sec,rn1);
char FILE_NAME[2048+1024];
 
// ###########################################################################
// ########################## Coupling Element _ Rotating Frame LA ###########
// ###########################################################################
complex double *la_coup     = calloc(2*N_b,sizeof(double));
        double *n_ph        = calloc(  N_b,sizeof(double));
complex double phon_fak[N_b];

INITIALIZE_PHONON_INTERACTION(la_coup,n_ph,phon_fak);

complex double phon_tab_n[N_b];
complex double phon_tab_m[N_b];
for (k=0;k<N_b;k++)
{
	phon_tab_n[k] = 1.;
	phon_tab_m[k] = 1.;
}

// ###########################################################################
// ########################## Coupling Element _ Rotating Frame LA ###########
// ###########################################################################

// ###########################################################################
// ########################## PATH_MAPPING FOR AUGMENTED DENSITY MATRIX ######
// ###########################################################################
// there are four options per time increment 0: 0 - l0,r0;  1 - l0,r1;  2 - l1,r0;  3 - l1,r1; 
int Nc=15; 
int pow4[Nc+1]; 
pow4[0]=1; for(i=1;i<Nc+1;i++) pow4[i] = 4*pow4[i-1];   
/*
4^ 1: 4        -- 4^ 2: 16        -- 4^ 3:     64 -- 4^ 4:     256 -- 4^ 5:    1024 -- 4^ 6:     4096 -- 
4^ 7: 16384    -- 4^ 8: 65536     -- 4^ 9: 262144 -- 4^10: 1048576 -- 4^11: 4194304 -- 4^12: 16777216 
4^13: 67108864 -- 4^14: 268435456 -- 4^15: 1073741824 */
int INIT_RHO = 0;
int MAX_K = 18200000;

int *path_table  = calloc( MAX_K*(Nc+1),sizeof(int));
CREATE_PATH_TABLE(path_table,MAX_K,Nc,pow4,INIT_RHO);
// with the path table I have an ordered series of all density matrices possible 
// first extract the number of possibitilies
int Pmax=0;
int INDEX[Nc];
i=1; m=1; INDEX[0]=0;
while(PATHS(i,0)>0) 
{
	if ( PATHS(i,0)==pow4[m] ) { INDEX[m]=i;m++; }
//	printf("PATHS(%i,0)=%i \n",i,PATHS(i,0));
	i++; 
} Pmax=i-1;
printf("PATHS up to %i with initial condition rho(0)=%i \n",Pmax,INIT_RHO);
for(i=0;i<m;i++) printf("MANIFOLDS P[%i]=%*i -- %*i \n",i,7,INDEX[i],7,PATHS(INDEX[i],0));
int MAN_MAX = m;
// ###########################################################################
// #################### END - PATH_MAPPING FOR AUGMENTED DENSITY MATRIX ######
// ###########################################################################

// ###########################################################################
// ########################## Dynamics #######################################
// ###########################################################################

// now calculate the electron-light density matrices 
complex double *adm = calloc( 2*MAX_K , sizeof(double) ); 

double OmR = creal(OmL);
complex double *Hsys        = calloc(2*20,sizeof(double));
// Define Bloch rotations
M_Om(0,0) = cos(OmR*delta_t)   ; M_Om(1,1) = M_Om(0,0); 
M_Om(0,1) = -I*sin(OmR*delta_t); M_Om(1,0) = M_Om(0,1); 
   
Md_Om(0,0)=conj(M_Om(0,0)); Md_Om(1,1)=conj(M_Om(1,1));
Md_Om(1,0)=conj(M_Om(0,1)); Md_Om(0,1)=conj(M_Om(1,0));


adm[0]=INIT_RHO; int D=1; int lpast,lnow,rpast,rnow; complex double ctemp,dtemp,etemp;
//printf("Re[adm]=%.10f -- Im[adm]=%.10f \n",creal(M_Om(1,0)*Md_Om(0,1)),cimag(M_Om(1,0)*Md_Om(0,1)));
//ctemp = M_Om(0,0)*Md_Om(0,0); 
//printf("Re[adm]=%.10f -- Im[adm]=%.10f \n",creal(ctemp*M_Om(1,0)*Md_Om(0,1)),cimag(ctemp*M_Om(1,0)*Md_Om(0,1)));
  
for (k=1;k<Pmax;k++)
{	
 if (PATHS(k,0)>pow4[D+1]) D++; // as many matrix multiplications as time steps given in powers of 4
// printf(" k=%i and Path=%i --> %i Matrix elements \n\n",k,PATHS(k,0),D); 
 adm[PATHS(k,0)] = 1.;
 for(i=1;i<D+1;i++)
 { 
  if ( PATHS(k,i)==0 ) { lpast=0;rpast=0; }
  if ( PATHS(k,i)==1 ) { lpast=0;rpast=1; }
  if ( PATHS(k,i)==2 ) { lpast=1;rpast=0; }
  if ( PATHS(k,i)==3 ) { lpast=1;rpast=1; }
 
  if ( PATHS(k,i+1)==0 ) { lnow=0;rnow=0; }
  if ( PATHS(k,i+1)==1 ) { lnow=0;rnow=1; }
  if ( PATHS(k,i+1)==2 ) { lnow=1;rnow=0; }
  if ( PATHS(k,i+1)==3 ) { lnow=1;rnow=1; }
  
  adm[PATHS(k,0)] *= M_Om(lnow,lpast) * Md_Om(rpast,rnow);
  //printf("k=%i:[%i]:(%i,%i) --> M_Om(%i,%i) * Md_OM(%i,%i) \n",k,PATHS(k,0),PATHS(k,i),PATHS(k,i+1),lnow,lpast,rpast,rnow); 
  //printf("k=%i:Re[adm]=%.10f -- Im[adm]=%.10f -(%i,%i,%i,%i, \n",k,creal(adm[k]),cimag(adm[k]),PATHS(k,1),PATHS(k,2),PATHS(k,3),PATHS(k,4));
 }
//printf("PATHS(%i,0)=%i \n",k,PATHS(k,0)); 
}

complex double rho11 = 0.;
rho11 = 0.; i=INDEX[1]; 
for(m=2;m<MAN_MAX;m++)
{	
rho11 = 0.;  
for(i=INDEX[m];i<INDEX[m+1];i++)
{ 	
 if ( PATHS(i,m)==3) rho11 += adm[PATHS(i,0)]; 
}	
printf("%i :ADM_11=%.10f -- EXACT=%.10f \n",m,creal(rho11),1.-cos(OmR*(m-1)*delta_t)*cos(OmR*(m-1)*delta_t) );
}

printf("\nADM_11=%.10f -- EXACT=%.10f \n",creal(adm[13]),1.-cos(OmR*(m-1)*delta_t)*cos(OmR*(m-1)*delta_t) );


i=0;m=1;
/*while( PATHS(m,0)>0 ) 
{
printf("STEP:%*i [%*i] (%i,%i,%i,%i \n",3,m,3,PATHS(m,0),PATHS(m,1),PATHS(m,2),PATHS(m,3),PATHS(m,4));
m++;
}*/
// time step [0]->[1] 00,01,02,03 
// time step [1]->[2] 0-00,0-01,0-02,0-03,0-10,0-11,0-12,0-13,0-20,0-21,0-22,0-23,0-30,0-31,0-32,0-33
complex double rho11_0=1.;
complex double rho11_1=0.;


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
//printf("Nc = %i \n",CUT_OFF);

complex double num_pol = 1.;
complex double ana_pol = 0.;

complex double IBM_arg = 0.;
int n,cut;	
int ehour = loctime -> tm_hour; int eminute = loctime -> tm_min; int esecond = loctime -> tm_sec;
int eyear = 1900+loctime -> tm_year; int emonth = loctime -> tm_mon; int eday  = loctime -> tm_mday;
double time_before,time_now; time_before =esecond + 60*eminute+60*60*ehour;

Snm = 0. + I*0.;

// initial condition


 /*
for(k=0; k<N; k++) // for increases the index after the loop is completed, k=0 in the first step
   {  
   t=delta_t*k;
  

 

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
free(path_table);
//free(adm);
free(n_ph);
return 0;
}

