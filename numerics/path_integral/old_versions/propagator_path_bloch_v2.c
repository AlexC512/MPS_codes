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

#define PATHS(m,n) ( *(path_table + (m)*(Nc+1) + (n) ) )
// rho(0,n) [      0--1*1048576-1, 1*4^10-1]
// rho(1,n) [1048576--2*1048576-1, 2*4^10-1]
// I don't care about the first time steps, where this density matrix is too large
#define rho(m,n)   ( *(adm + (m)*1048576 + (n) )) 

/*
4^ 1: 4        -- 4^ 2: 16        -- 4^ 3:     64 -- 4^ 4:     256 -- 4^ 5:    1024 -- 4^ 6:     4096 -- 
4^ 7: 16384    -- 4^ 8: 65536     -- 4^ 9: 262144 -- 4^10: 1048576 -- 4^11: 4194304 -- 4^12: 16777216 
4^13: 67108864 -- 4^14: 268435456 -- 4^15: 1073741824 */

		
// #################### Parameter-Uebernahme Ende ############################
//////////////////////////////////////////////////////////////////////////////
// ################### normalized transition coupling element phonon #########

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
// input integer k and array of length MAX_K and 4powers, output array (0,0,1, .. dependent on k
void CREATE_PATH_VECTOR(int *path_vector,int INT_PATH, int Nc, int *pow4)
{

int tPATH[Nc]; // temp vector
int i,k_map;
// input nPATH, I start with 
k_map = INT_PATH;
i = Nc; 	
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
for (i=0;i<Nc;i++) path_vector[i]=tPATH[i];

return;
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
// I need an augmented density matrix for timestep N, with 4^N+1 elements, I will start with 
// an adm of ten time steps ...  time step zero adm[0-3] ... time step one adm[0-15]  ...   
// how do I map rho(0,0-3), rho(1,0-15), rho(2,0-63) .... rho(10,4^11) 
// rho(1) hat sechzehn einträge, ich muss jeder zahl 0 0000, 1, 0001, etc einen path zu ordnen,
// dafür brauche ich funktion die aus k --> (0,0,...1, ) einen vektor mit einträgen macht
// dann lasse ich loop von 0-4^k laufen, und ordne jedem integer ein path zu, und dann schreibe 
// ich direkt den bloch propagator! ACHTUNG, ich kann nicht länger werden als 4^15, also muss ich 
// ein Tensorprodukt definieren für längere Memoryketten 

int Nc=10; 
int pow4[Nc+1]; 
pow4[0]=1; for(i=1;i<Nc+1;i++) pow4[i] = 4*pow4[i-1];   
/*
4^ 1: 4        -- 4^ 2: 16        -- 4^ 3:     64 -- 4^ 4:     256 -- 4^ 5:    1024 -- 4^ 6:     4096 -- 
4^ 7: 16384    -- 4^ 8: 65536     -- 4^ 9: 262144 -- 4^10: 1048576 -- 4^11: 4194304 -- 4^12: 16777216 
4^13: 67108864 -- 4^14: 268435456 -- 4^15: 1073741824 */

// mapping an integer from 0 - 4^10-1 to sequence of m_i=(0,1,2,3) (m0,m1,m2,m3,m4,m5,m6,m7,m8,m9) 
// test run
int PATH[Nc];
CREATE_PATH_VECTOR(PATH,1048576-1,Nc,pow4);

printf("(%i",PATH[0]);  for (i=1;i<Nc;i++) printf(",%i",PATH[i]); printf(")\n"); 

// ###########################################################################
// #################### END - PATH_MAPPING FOR AUGMENTED DENSITY MATRIX ######
// ###########################################################################

// ###########################################################################
// ########################## Dynamics #######################################
// ###########################################################################

// the augmented density matrix for a memory length of ten consumes 10x4^10 complex numbers  
complex double *adm = calloc( 2*pow4[10]*10 , sizeof(double) ); 

double OmR = creal(OmL);
complex double *Hsys        = calloc(2*20,sizeof(double));
// Define Bloch rotations
M_Om(0,0) = cos(OmR*delta_t)   ; M_Om(1,1) = M_Om(0,0); 
M_Om(0,1) = -I*sin(OmR*delta_t); M_Om(1,0) = M_Om(0,1); 
   
Md_Om(0,0)=conj(M_Om(0,0)); Md_Om(1,1)=Md_Om(1,1);
Md_Om(1,0)=conj(M_Om(0,1)); Md_Om(0,1)=Md_Om(1,0);


int D=1; int lpast,lnow,rpast,rnow; 
// -------------------
rho(0,0)=1; //ground state populated
rho(0,1)=0; //ground state populated
rho(0,2)=0; //ground state populated
rho(0,3)=0; //ground state populated
// --------------------
for (k=0;k<pow4[2];k++)
{	
 D=1;
 rho(1,k) = 1.;
 CREATE_PATH_VECTOR(PATH,k,Nc,pow4);
 for(i=0;i<D;i++)
 { 
  if ( PATH[i]==0 ) { lpast=0;rpast=0; }
  if ( PATH[i]==1 ) { lpast=0;rpast=1; }
  if ( PATH[i]==2 ) { lpast=1;rpast=0; }
  if ( PATH[i]==3 ) { lpast=1;rpast=1; }
 
  if ( PATH[i+1]==0 ) { lnow=0;rnow=0; }
  if ( PATH[i+1]==1 ) { lnow=0;rnow=1; }
  if ( PATH[i+1]==2 ) { lnow=1;rnow=0; }
  if ( PATH[i+1]==3 ) { lnow=1;rnow=1; }
  // implement the initial density matrix by hand
  if ( (i==0) && (lpast==0) && (rpast==0) ) { if ( (lnow==0) && (rnow==0) ) rho(1,k) *=     M_Om(0,0)*M_Om(0,0);
  											  if ( (lnow==1) && (rnow==1) ) rho(1,k) *= (1.-M_Om(0,0)*M_Om(0,0));		
                                              if ( (lnow==0) && (rnow==1) ) rho(1,k) *= M_Om(0,0)*Md_Om(0,1);  
                                              if ( (lnow==1) && (rnow==0) ) rho(1,k) *= Md_Om(0,0)*M_Om(1,0);   
											}
  else { if (i==0) rho(1,k) *= 0;
         else      rho(1,k) *= M_Om(lnow,lpast) * Md_Om(rpast,rnow);
       }	 
  //printf("k=%i:(%i,%i) --> M_Om(%i,%i) * Md_OM(%i,%i) \n",k,PATH[i],PATH[i+1],lnow,lpast,rpast,rnow); 
 //printf("k=%i:Re[adm]=%.10f -- Im[adm]=%.10f -(%i,%i,%i, \n",k,creal(rho(1,k)),cimag(rho(1,k)),PATH[0],PATH[1],PATH[2]);
 }
}
// expectation value take trace
complex double rho11=0.;
for (k=0;k<pow4[2];k++)
{	
 CREATE_PATH_VECTOR(PATH,k,Nc,pow4);
 
 printf("k=%i:Re[adm]=%.10f -- Im[adm]=%.10f -(%i,%i,%i, \n",k,creal(rho(1,k)),cimag(rho(1,k)),PATH[0],PATH[1],PATH[2]); 
 if (PATH[1]==3) rho11 += rho(1,k);
}
printf("\nADM_11=%.10f -- EXACT=%.10f \n",creal(rho11),1.-cos(OmR*delta_t)*cos(OmR*delta_t) );
// --------------------
for (k=0;k<pow4[3];k++)
{	
 D=2;
 rho(2,k) = 1.;
 CREATE_PATH_VECTOR(PATH,k,Nc,pow4);
 for(i=0;i<D;i++)
 { 
  if ( PATH[i]==0 ) { lpast=0;rpast=0; }
  if ( PATH[i]==1 ) { lpast=0;rpast=1; }
  if ( PATH[i]==2 ) { lpast=1;rpast=0; }
  if ( PATH[i]==3 ) { lpast=1;rpast=1; }
 
  if ( PATH[i+1]==0 ) { lnow=0;rnow=0; }
  if ( PATH[i+1]==1 ) { lnow=0;rnow=1; }
  if ( PATH[i+1]==2 ) { lnow=1;rnow=0; }
  if ( PATH[i+1]==3 ) { lnow=1;rnow=1; }
  // implement the initial density matrix by hand
  if ( (i==0) && (lpast==0) && (rpast==0) ) { if ( (lnow==0) && (rnow==0) ) rho(2,k) *=     M_Om(0,0)*M_Om(0,0);
  											  if ( (lnow==1) && (rnow==1) ) rho(2,k) *= (1.-M_Om(0,0)*M_Om(0,0));		
                                              if ( (lnow==0) && (rnow==1) ) rho(2,k) *= M_Om(0,0)*Md_Om(0,1);  
                                              if ( (lnow==1) && (rnow==0) ) rho(2,k) *= Md_Om(0,0)*M_Om(1,0);   
											if (k==52) printf("yes  %.10f\n",cimag(rho(2,k)));
											}
  else { if (i==0) {rho(2,k) *= 0;  if (k==52) printf("%i yes  %.10f\n",i,creal(rho(2,k)));
							}
         else     { rho(2,k) *= M_Om(lnow,lpast) ;
		            rho(2,k) *= Md_Om(rpast,rnow);
		 
		 }
		 
		 			if (k==52) printf("%i yes  %.10f\n",i,cimag(rho(2,k)));
								
       }	 
  printf("k=%i:(%i,%i) --> M_Om(%i,%i) * Md_OM(%i,%i) \n",k,PATH[i],PATH[i+1],lnow,lpast,rpast,rnow); 
 printf("k=%i:Re[adm]=%.10f -- Im[adm]=%.10f -(%i,%i,%i, \n",k,creal(rho(2,k)),cimag(rho(2,k)),PATH[0],PATH[1],PATH[2]);
 }
}
// expectation value take trace
rho11=0.;
for (k=0;k<pow4[3];k++)
{	
 CREATE_PATH_VECTOR(PATH,k,Nc,pow4);
 printf("k=%i:Re[adm]=%.10f -- Im[adm]=%.10f -(%i,%i,%i, \n",k,creal(rho(2,k)),cimag(rho(2,k)),PATH[0],PATH[1],PATH[2]); 
 
 if (PATH[2]==3) rho11 += rho(2,k);
}
printf("\nADM_11=%.10f -- EXACT=%.10f \n",1.*creal(rho11),1.-cos(OmR*2*delta_t)*cos(OmR*2*delta_t) );






//free(path_table);
free(adm);
free(la_coup);
free(n_ph);
return 0;
}

