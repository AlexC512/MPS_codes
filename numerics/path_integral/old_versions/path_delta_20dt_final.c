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

 #define TIME_DYNAMICS_TO_NS 0.10000
 #define T_STEPS 250
 #define CUT_OFF_DENSITY ( 0.00000001 * 10.01  )
 #define N_DELTA_T ( 20 )
 int N=T_STEPS;  
 double delta_t = ((double)TIME_DYNAMICS_TO_NS) * 1000. * 1000. /((double)T_STEPS); 
 #define MEMORY_CUT (20.0*1000.0)
 // ############# ELECTRONIC PARAMETERS (start) ########################################################
 #define GAM_R (0.000000633 *0.) // (fs)^(-1) ... T1 time 750ps=750 000fs, GAM_R = 1/2T1 = 0.000000666
 // (fs)^(-1) ... T1 time 600ps=600 000fs, GAM_R = 1/2T1 = 0.000000833
 // (fs)^(-1) ... T1 time 537ps=537 000fs, GAM_R = 1/2T1 = 0.000000933
 #define AW_VC (0.0)
 #define AW_CC (0.0)
 #define AW_VV (1.0)
 #define AW_CV (0.0)
 
 #define OmL      ( 0.0000784 * 0.0  ) 
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
    #define rho_mat (33514170.0) // eVfs^2 /nm^5
    // kg = ... nutze Energie=Kraft x Weg , [kg] = 6.24151*10^3 ev fs^2/nm^5  
    // eV =1.6021766208(98)×10^19 J 
    // rho 5370 kg/m^3 --> 5370*6241,51 eV fs^2/nm^5 = 33514170 eV fs^2/nm^5
    #define D_c (-14.6) // eV
    #define D_v (-4.8)  // eV
 
// ############# PHONONIC PARAMETERS (end) ############################################################
// ############# RK4 PARAMETERS (end) ###################################################################

#define M_Om(m,n)    ( *(Hsys      + 0 + 2*(m) + (n) )) 
#define Md_Om(m,n)   ( *(Hsys      + 4 + 2*(m) + (n) )) 

// rho(0,n) [      0--1*1048576-1, 1*4^10-1]
// rho(1,n) [1048576--2*1048576-1, 2*4^10-1]
// I don't care about the first time steps, where this density matrix is too large
#define rs_index(n,p)     ( *(REL_STATES_INDEX   + (n)*2 + (p) ) ) 
#define rs_2_index(n,p)   ( *(REL_STATES_INDEX_2 + (n)*2 + (p) ) ) 


#define rho_past(m,n)   ( *(adm_past + (m)*4 + (n) )) 

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
// input integer k and array of length MAX_K and 4powers, output array (0,0,1, .. dependent on k
void CREATE_DOUBLE_PATH_VECTOR(int *path_vector,int INT_PATH1,int INT_PATH2, int MAX_POWER, int *pow4)
{

int tPATH[MAX_POWER]; // temp vector
int i,k_map;
// input nPATH, I start with 
k_map = INT_PATH1;
i = MAX_POWER; 	
	 while( (k_map>0) && (i>=0) )  
	 { 
        if ( pow4[i] > k_map ) tPATH[i]=0; 
		else 			
		{
	     tPATH[i] = k_map/pow4[i] ;
		 k_map   -= pow4[i]*tPATH[i];
     	}
	 i--; 	
     }
	 // here I fill in the zeros 
	 while ( i>=0 ) 
	 { 
       tPATH[i]=0; 
	   i--;
	 }

//if max power 11, than 0-10 is the first part 
for (i=0;i<MAX_POWER;i++)  { path_vector[i]=tPATH[i]; tPATH[i]=0;  }

// input nPATH, I start with 
k_map = INT_PATH2;
i = MAX_POWER; 	
	 while( (k_map>0) && (i>=0) )  
	 { 
        if ( pow4[i] > k_map ) tPATH[i]=0; 
		else 			
		{
	     tPATH[i] = k_map/pow4[i] ;
		 k_map   -= pow4[i]*tPATH[i];
     	}
	 i--; 	
     }
	 // here I fill in the zeros 
	 while ( i>=0 ) 
	 { 
       tPATH[i]=0; 
	   i--;
	 }

//printf("(%i",tPATH[0]);  for (i=1;i<Nc;i++) printf(",%i",tPATH[i]); printf(")\n"); 
for (i=0;i<MAX_POWER;i++) path_vector[i+MAX_POWER]=tPATH[i];

return;
}
// #######################################################
// #######################################################
// #######################################################

void PRINT_PATH(int *PATH, int Nc)
{
int i; 
printf("(%*i",2,PATH[0]);  for (i=1;i<Nc+1;i++) printf(",%*i",2,PATH[i]); printf(")\n"); 
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
// I need an augmented density matrix for timestep N, with 4^N+1 elements, I will start with 
// an adm of ten time steps ...  time step zero adm[0-3] ... time step one adm[0-15]  ...   
// how do I map rho(0,0-3), rho(1,0-15), rho(2,0-63) .... rho(10,4^11) 
// rho(1) hat sechzehn einträge, ich muss jeder zahl 0 0000, 1, 0001, etc einen path zu ordnen,
// dafür brauche ich funktion die aus k --> (0,0,...1, ) einen vektor mit einträgen macht
// dann lasse ich loop von 0-4^k laufen, und ordne jedem integer ein path zu, und dann schreibe 
// ich direkt den bloch propagator! ACHTUNG, ich kann nicht länger werden als 4^15, also muss ich 
// ein Tensorprodukt definieren für längere Memoryketten 

int Nc=21; 
int MAX_POWER=11;
int pow4[MAX_POWER]; // we don't go higher due to number length issues 
pow4[0]=1; for(i=1;i<MAX_POWER+1;i++) pow4[i] = 4*pow4[i-1];   
PRINT_PATH(pow4,MAX_POWER);
int PATH[Nc]; 
int PATH_INDEX[3];
// actually length of path 20, but we need 21 one steps 
// checking whether mapping works
printf("(%*i",2,0);  for (i=1;i<Nc+1;i++) printf(",%*i",2,i); printf(")\n"); 
CREATE_DOUBLE_PATH_VECTOR(PATH,pow4[11]-1,0,MAX_POWER,pow4); 
PRINT_PATH(PATH,Nc);
CREATE_DOUBLE_PATH_VECTOR(PATH,pow4[11]-1,pow4[11]-1,MAX_POWER,pow4); 
PRINT_PATH(PATH,Nc);
CREATE_DOUBLE_PATH_VECTOR(PATH,0,1,MAX_POWER,pow4); 
PRINT_PATH(PATH,Nc);
printf("PATH[0] -- PATH[10]...PATH[11] -- PATH[21] \n");
printf("POW4[0] -- POW4[10]...POW4[ 0] -- POW4[10] \n");
// ###########################################################################
// #################### END - PATH_MAPPING FOR AUGMENTED DENSITY MATRIX ######
// ###########################################################################

complex double *adm_past = calloc( 2*4*Nc        , sizeof(double) ); 

// ############ Define Bloch rotations #######################
double OmR = creal(OmL);
complex double *Hsys        = calloc(2*20,sizeof(double));
// delta-pulse conditions !!!!!!!
M_Om(0,0) = 1.   ; M_Om(1,1) = M_Om(0,0); 
M_Om(0,1) = -I*0.5; M_Om(1,0) = M_Om(0,1); 
Md_Om(0,0)=M_Om(0,0); Md_Om(1,1)=M_Om(1,1);
Md_Om(1,0)=conj(M_Om(0,1)); Md_Om(0,1)=Md_Om(1,0);
// ############ End Define Bloch rotations #######################

int D=1; int lpast,lnow,rpast,rnow,p,l; 
complex double rho11=0.;
complex double rho10=0.;
complex double rho01=0.;
complex double rho00=0.;

// -------------------
for(p=0;p<Nc;p++)
{
rho_past(p,0)=1; //ground state populated
rho_past(p,1)=1; //ground state populated
rho_past(p,2)=0; //ground state populated
rho_past(p,3)=0; //ground state populated
}
// --------------------
complex double init_rho;
complex double test=1.;
int NEGLIGIBLE =0;
int NON_NEGLIGIBLE = 0;
double OmRp;
int ind0,ind1,L_LIMIT,K_LIMIT;
int MAX_NN;

//int REL_STATES[65536]; no safe of adm yet
int *REL_STATES_INDEX   = calloc( 2*(16777216+1), sizeof(int));
int *REL_STATES_INDEX_2 = calloc( 2*(16777216+1), sizeof(int));

// ###########################################################################
// ###########################################################################
// ########################## DYNAMICS -- 0 --> 11 DELTA_T ###################
// ###########################################################################
// ###########################################################################


// I start with calculating the first Nc time steps (the memory depth)
// Afterwards, I need to repeat this Nc time steps to reach the next step but 
// with an incremental initial condition -- this is the basic problem of the 
// non Markovian problem 
complex double init_00=1.;   
complex double init_01=0.;   
complex double init_10=0.;   
complex double init_11=0.;   

complex double temp_00=0.;   
complex double temp_01=0.;   
complex double temp_10=0.;   
complex double temp_11=0.;   

// PATH[0] INITIAL STATE, p counts time steps 
for(p=1;p<5;p++){	
rho11=0.;rho10=0.;rho01=0.;rho00=0.;
if (p<=10) { L_LIMIT=1;          K_LIMIT = pow4[p+1]; } 
else       { L_LIMIT=pow4[p-10]; K_LIMIT = pow4[ 11]; }

for (l=0;l<L_LIMIT;l++){
for (k=0;k<K_LIMIT;k++)
{	
 test = 1.;
 CREATE_DOUBLE_PATH_VECTOR(PATH,k,l,MAX_POWER,pow4);
 //if ((l>0) && (k%100000==0) ) PRINT_PATH(PATH,Nc);
 for(i=0;i<p;i++)
 {

  if ( PATH[i]==0 ) { lpast=0;rpast=0; }
  if ( PATH[i]==1 ) { lpast=0;rpast=1; }
  if ( PATH[i]==2 ) { lpast=1;rpast=0; }
  if ( PATH[i]==3 ) { lpast=1;rpast=1; }
 
  if ( PATH[i+1]==0 ) { lnow=0;rnow=0; }
  if ( PATH[i+1]==1 ) { lnow=0;rnow=1; }
  if ( PATH[i+1]==2 ) { lnow=1;rnow=0; }
  if ( PATH[i+1]==3 ) { lnow=1;rnow=1; }
  // include the initial state of the adm
    if (i==0)   
       { 
        if ((lpast==1)&&(rpast==1)) init_rho  = init_11;
        if ((lpast==0)&&(rpast==1)) init_rho  = init_01;
        if ((lpast==1)&&(rpast==0)) init_rho  = init_10;
        if ((lpast==0)&&(rpast==0)) init_rho  = init_00;
        // for the delta-pulse only at i=0 there is an off-diagonal contribution
        test *= M_Om(lnow,lpast) * init_rho * Md_Om(rpast,rnow); 
	   }
	// if i>0 no light coupling, just identity   
    else test *= 1.;
 }
 
 if ( ( cabs(test) > CUT_OFF_DENSITY ) )   
    {
        if (PATH[p]==0) rho00 += test;
        if (PATH[p]==1) rho01 += test;
        if (PATH[p]==2) rho10 += test;
		if (PATH[p]==3) rho11 += test;
		
        if (p==4) 
		{ 
         ind0=0;ind1=0;
		 for(m=0;m<MAX_POWER;m++) { 
		                    ind0  +=  pow4[m]*PATH[m+ 0];
		                    ind1  +=  pow4[m]*PATH[m+MAX_POWER]; 
		                   }
         rs_2_index(NON_NEGLIGIBLE,0) = ind0;
		 rs_2_index(NON_NEGLIGIBLE,1) = ind1;		
         NON_NEGLIGIBLE++;
		} 

     }
} // k loop
} // l loop

if (p==1) 
{
	temp_00 = rho00;
    temp_01 = rho01;
	temp_10 = rho10;
    temp_11 = rho11;
}	

OmRp = p*delta_t*OmR;
printf("[%i]:t=%.5fps: ADM_00=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",p,OmRp/OmR*0.001,creal(rho00),cos(OmRp)*cos(OmRp),creal(rho00)-cos(OmRp)*cos(OmRp) );
printf("[%i]:t=%.5fps: ADM_01=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",p,OmRp/OmR*0.001,cimag(rho01),sin(OmRp)*cos(OmRp),cimag(rho01)-sin(OmRp)*cos(OmRp) );
printf("[%i]:t=%.5fps: ADM_10=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",p,OmRp/OmR*0.001,cimag(rho10),-sin(OmRp)*cos(OmRp),cimag(rho10)+sin(OmRp)*cos(OmRp) );
printf("[%i]:t=%.5fps: ADM_11=%.10f -- EXACT=%.10f -- DIFF=%.10f \n\n",p,OmRp/OmR*0.001,creal(rho11),sin(OmRp)*sin(OmRp),creal(rho11)-sin(OmRp)*sin(OmRp) );
//PRINT_PATH(PATH,Nc);
}

printf("RELEVANT STATES %i \n",NON_NEGLIGIBLE);
/*
// ###########################################################################
// ###########################################################################
// ########################## DYNAMICS -- 12 --> 21 DELTA_T ##################
// ###########################################################################
// ###########################################################################
    
for(p=5;p<22;p++)
{

MAX_NN = NON_NEGLIGIBLE; 	
for(k=0;k<MAX_NN;k++) 
   { 
    ind0 = rs_2_index(k,0);
    ind1 = rs_2_index(k,1);
    rs_index(k,0) = ind0; 
	rs_index(k,1) = ind1; 
	REL_STATES_INDEX_2[k]=0;
   } 
NON_NEGLIGIBLE=0;	

rho11=0.;rho10=0.;rho01=0.;rho00=0.;
for(l=0;l<4;l++)
{ 
for(k=0;k<MAX_NN;k++)
{	

 test = 1.;
 CREATE_DOUBLE_PATH_VECTOR(PATH,rs_index(k,0),rs_index(k,1),MAX_POWER,pow4);
 //PRINT_PATH(PATH,Nc);
 PATH[p]=l; 
 
 for(i=0;i<p;i++)
 {
  if ( PATH[i]==0 ) { lpast=0;rpast=0; }
  if ( PATH[i]==1 ) { lpast=0;rpast=1; }
  if ( PATH[i]==2 ) { lpast=1;rpast=0; }
  if ( PATH[i]==3 ) { lpast=1;rpast=1; }
 
  if ( PATH[i+1]==0 ) { lnow=0;rnow=0; }
  if ( PATH[i+1]==1 ) { lnow=0;rnow=1; }
  if ( PATH[i+1]==2 ) { lnow=1;rnow=0; }
  if ( PATH[i+1]==3 ) { lnow=1;rnow=1; }
  // include the initial state of the adm
    if (i==0)   
       { 
        if ((lpast==1)&&(rpast==1)) init_rho  = init_11;
        if ((lpast==0)&&(rpast==1)) init_rho  = init_01;
        if ((lpast==1)&&(rpast==0)) init_rho  = init_10;
        if ((lpast==0)&&(rpast==0)) init_rho  = init_00;
        // initial delta pulse ... valid only at i==0
        test *= 1.;  
       }
    else test *= 1.;
 }
 
 if ( ( cabs(test) > CUT_OFF_DENSITY ) )
    {
        if (PATH[p]==0) rho00 += test;
        if (PATH[p]==1) rho01 += test;
        if (PATH[p]==2) rho10 += test;
        if (PATH[p]==3) rho11 += test;
       
         ind0=0;ind1=0;
		 for(m=0;m<MAX_POWER;m++) { 
		                    ind0  +=  pow4[m]*PATH[m+ 0];
		                    ind1  +=  pow4[m]*PATH[m+MAX_POWER]; 
		                   }
         rs_2_index(NON_NEGLIGIBLE,0) = ind0;
		 rs_2_index(NON_NEGLIGIBLE,1) = ind1;		
         NON_NEGLIGIBLE++;
    }
 }    
} 
// Benchmark with exact result 
OmRp = p*delta_t*OmR;
printf("p=%i:t=%.5fps: ADM_00=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",p,OmRp/OmR*0.001,creal(rho00),cos(OmRp)*cos(OmRp),creal(rho00)-cos(OmRp)*cos(OmRp) );
printf("p=%i:t=%.5fps: ADM_01=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",p,OmRp/OmR*0.001,cimag(rho01),sin(OmRp)*cos(OmRp),cimag(rho01)-sin(OmRp)*cos(OmRp) );
printf("p=%i:t=%.5fps: ADM_10=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",p,OmRp/OmR*0.001,cimag(rho10),-sin(OmRp)*cos(OmRp),cimag(rho10)+sin(OmRp)*cos(OmRp) );
printf("p=%i:t=%.5fps: ADM_11=%.10f -- EXACT=%.10f -- DIFF=%.10f \n\n",p,OmRp/OmR*0.001,creal(rho11),sin(OmRp)*sin(OmRp),creal(rho11)-sin(OmRp)*sin(OmRp) );
printf("---- RELEVANT STATES: %i --- \n",NON_NEGLIGIBLE);
//PRINT_PATH(PATH,Nc);
}

printf("RELEVANT STATES %i \n",NON_NEGLIGIBLE);

// ###########################################################################
// ###########################################################################
// ###########################################################################

// After the first 20 time steps, we have to repeat them on and on for the various next time step, that is the memory problem
// no way to use the relevant states for a different initial conditions -- I guess. 
// the first steps are the fastest, as the memory is the shortest

printf("NOW WE HAVE TO REPEAT FROM STEP 1 TO GO TO STEP 22 \n");

// ###########################################################################
// ###########################################################################
// ########################## DYNAMICS -- from 21 DELTA_T on #################
// ###########################################################################
// ###########################################################################

int Dt;
for(Dt=1;Dt<N_DELTA_T;Dt++)
{

for (m=0;m<NON_NEGLIGIBLE;m++)
{
	REL_STATES_INDEX[m]=0;
	REL_STATES_INDEX_2[m]=0;
}
NON_NEGLIGIBLE =0;

init_00=temp_00;   
init_01=temp_01;   
init_10=temp_10;   
init_11=temp_11;   

for(p=1;p<5;p++){	
rho11=0.;rho10=0.;rho01=0.;rho00=0.;
if (p<=10) { L_LIMIT=1;          K_LIMIT = pow4[p+1]; } 
else       { L_LIMIT=pow4[p-10]; K_LIMIT = pow4[ 11]; }

for (l=0;l<L_LIMIT;l++){
for (k=0;k<K_LIMIT;k++)
{	
 test = 1.;
 CREATE_DOUBLE_PATH_VECTOR(PATH,k,l,MAX_POWER,pow4);

 for(i=0;i<p;i++)
 {

  if ( PATH[i]==0 ) { lpast=0;rpast=0; }
  if ( PATH[i]==1 ) { lpast=0;rpast=1; }
  if ( PATH[i]==2 ) { lpast=1;rpast=0; }
  if ( PATH[i]==3 ) { lpast=1;rpast=1; }
 
  if ( PATH[i+1]==0 ) { lnow=0;rnow=0; }
  if ( PATH[i+1]==1 ) { lnow=0;rnow=1; }
  if ( PATH[i+1]==2 ) { lnow=1;rnow=0; }
  if ( PATH[i+1]==3 ) { lnow=1;rnow=1; }
  // include the initial state of the adm
    if (i==0)   
       { 
        if ((lpast==1)&&(rpast==1)) init_rho  = init_11;
        if ((lpast==0)&&(rpast==1)) init_rho  = init_01;
        if ((lpast==1)&&(rpast==0)) init_rho  = init_10;
        if ((lpast==0)&&(rpast==0)) init_rho  = init_00;
        // the initial delta_pulse is now not active anymore
        test *= init_rho;  
       }
    else test *= 1.;
 }
 
 if ( ( cabs(test) > CUT_OFF_DENSITY ) )   
    {
        if (PATH[p]==0) rho00 += test;
        if (PATH[p]==1) rho01 += test;
        if (PATH[p]==2) rho10 += test;
		if (PATH[p]==3) rho11 += test;
		
        if (p==4) 
		{ 
         ind0=0;ind1=0;
		 for(m=0;m<MAX_POWER;m++) { 
		                    ind0  +=  pow4[m]*PATH[m+ 0];
		                    ind1  +=  pow4[m]*PATH[m+MAX_POWER]; 
		                   }
         rs_2_index(NON_NEGLIGIBLE,0) = ind0;
		 rs_2_index(NON_NEGLIGIBLE,1) = ind1;		
         NON_NEGLIGIBLE++;
		} 

     }
} // k loop
} // l loop

if (p==1)
{
temp_00 = rho00 ;
temp_01 = rho01 ;
temp_10 = rho10 ;
temp_11 = rho11 ;
}
}
    
for(p=5;p<22;p++)
{

MAX_NN = NON_NEGLIGIBLE; 	
for(k=0;k<MAX_NN;k++) 
   { 
    ind0 = rs_2_index(k,0);
    ind1 = rs_2_index(k,1);
    rs_index(k,0) = ind0; 
	rs_index(k,1) = ind1; 
	REL_STATES_INDEX_2[k]=0;
   } 
NON_NEGLIGIBLE=0;	

rho11=0.;rho10=0.;rho01=0.;rho00=0.;
for(l=0;l<4;l++)
{ 
for(k=0;k<MAX_NN;k++)
{	

 test = 1.;
 CREATE_DOUBLE_PATH_VECTOR(PATH,rs_index(k,0),rs_index(k,1),MAX_POWER,pow4);

 PATH[p]=l; 
 
 for(i=0;i<p;i++)
 {
  if ( PATH[i]==0 ) { lpast=0;rpast=0; }
  if ( PATH[i]==1 ) { lpast=0;rpast=1; }
  if ( PATH[i]==2 ) { lpast=1;rpast=0; }
  if ( PATH[i]==3 ) { lpast=1;rpast=1; }
 
  if ( PATH[i+1]==0 ) { lnow=0;rnow=0; }
  if ( PATH[i+1]==1 ) { lnow=0;rnow=1; }
  if ( PATH[i+1]==2 ) { lnow=1;rnow=0; }
  if ( PATH[i+1]==3 ) { lnow=1;rnow=1; }
  // include the initial state of the adm
    if (i==0)   
       { 
        if ((lpast==1)&&(rpast==1)) init_rho  = init_11;
        if ((lpast==0)&&(rpast==1)) init_rho  = init_01;
        if ((lpast==1)&&(rpast==0)) init_rho  = init_10;
        if ((lpast==0)&&(rpast==0)) init_rho  = init_00;
        // delta_pulse not valid anymore as i==0 is now delta_t
        test *= init_rho;  
       }
    else test *= 1.;
 }
 
 if ( ( cabs(test) > CUT_OFF_DENSITY ) )
    {
        if (PATH[p]==0) rho00 += test;
        if (PATH[p]==1) rho01 += test;
        if (PATH[p]==2) rho10 += test;
        if (PATH[p]==3) rho11 += test;
       
         ind0=0;ind1=0;
		 for(m=0;m<MAX_POWER;m++) { 
		                    ind0  +=  pow4[m]*PATH[m+ 0];
		                    ind1  +=  pow4[m]*PATH[m+MAX_POWER]; 
		                   }
         rs_2_index(NON_NEGLIGIBLE,0) = ind0;
		 rs_2_index(NON_NEGLIGIBLE,1) = ind1;		
         NON_NEGLIGIBLE++;
    }
 }    
} 
}

// Benchmark with exact result
 
OmRp = (Dt+21)*delta_t*OmR;
printf("p=%i:t=%.5fps: ADM_00=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",Dt+21,OmRp/OmR*0.001,creal(rho00),cos(OmRp)*cos(OmRp),creal(rho00 )-cos(OmRp)*cos(OmRp) );
printf("p=%i:t=%.5fps: ADM_01=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",Dt+21,OmRp/OmR*0.001,cimag(rho01),sin(OmRp)*cos(OmRp),cimag(rho01 )-sin(OmRp)*cos(OmRp) );
printf("p=%i:t=%.5fps: ADM_10=%.10f -- EXACT=%.10f -- DIFF=%.10f \n",Dt+21,OmRp/OmR*0.001,cimag(rho10 ),-sin(OmRp)*cos(OmRp),cimag(rho10)+sin(OmRp)*cos(OmRp) );
printf("p=%i:t=%.5fps: ADM_11=%.10f -- EXACT=%.10f -- DIFF=%.10f \n\n",Dt+21,OmRp/OmR*0.001,creal(rho11),sin(OmRp)*sin(OmRp),creal(rho11)-sin(OmRp)*sin(OmRp) );
printf("---- RELEVANT STATES: %i --- \n",NON_NEGLIGIBLE);

// ###########################################################################
// ###########################################################################
// ###########################################################################

} // Dt loop
*/

free(REL_STATES_INDEX);
free(REL_STATES_INDEX_2);
free(la_coup);
free(n_ph);
return 0;
}
