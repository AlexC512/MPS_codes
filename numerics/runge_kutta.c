#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <time.h>
#define hbar ( 0.658212196 )     // eVfs
#define PI M_PI
// ################ ZEITRAUM #######################
#define TIME_DYNAMICS_TO ( 0.01000 )
#define T_STEPS 50000
int N     = T_STEPS;
int STEP_IS_SAVED= 1; // jeder ... Schritt wird in ein File geschrieben!!
int TIME_OUTPUT  =2000;
double delta_t = ((double)TIME_DYNAMICS_TO) * 1000. * 1000. /( (double)T_STEPS);

#define GAM ( 0.23   * 1. ) //( Mc/50)
#define gam_p ( 0.500   * 1. ) //( Mc/50)
#define OM ( 2.10003   * 1. ) //( Mc/50)

// #########################################################
// ############# RK4 #######################################
// #########################################################

#define CALCULATE_TEMP(a,b,cb,dest,i) ({for(i=0;i<N_DGL;i++) (dest)[i]=(a)[i]+(cb)*(b)[i];})
#define ADD_VECTOR(dest,in,in_fak,i) ({for(i=0;i<N_DGL;i++) dest[i]=dest[i]+(in_fak)* (in)[i];})

#define SECHSTEL 0.1666666666666666666666667
#define DRITTEL  0.3333333333333333333333333
#define HALFTE   0.5
// Zur Definition der Zustaende
#define N_l (4) // Zahl der Level
#define rho(a,b)      ( *(derivates     + (a)*N_l + (b)      ))
#define rho_out(a,b)  ( *(derivates_out + (a)*N_l + (b)      ))
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
   OUTPUT[99] += cabs(derivates[i]-deriv_stst[i])*delta_t;
   deriv_stst[i] = derivates[i];
 }
 return;
}
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

double pulse = OM*exp(-t*t*100.0);

for (a=0;a<100;a++) OUTPUT[a]=0.;
 
rho_out(0,0) =  +2.*GAM           *rho(1,1)+I*OM*( rho(0,1)-rho(1,0) );
rho_out(1,0) = - 2.*(GAM/2.+gam_p)*rho(1,0)+I*OM*( rho(1,1)-rho(0,0) );
rho_out(0,1) = - 2.*(GAM/2.+gam_p)*rho(0,1)-I*OM*( rho(1,1)-rho(0,0) );
rho_out(1,1) = - 2.*GAM           *rho(1,1)-I*OM*( rho(0,1)-rho(1,0) );

  return ;
}

// ###############################################################################
// ###############################################################################

void RUNGE_KUTTA
( complex double *derivates,
  complex double *k1,
  complex double *k2,
  complex double *temp,
  double *OUTPUT,
  double t,
  int N_DGL
)
{
 int i;
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

 return;
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

FILE *f_spectrum;
f_spectrum=fopen("dummy_2.dat","w");
// ###############################################################################
// ############################## INITIAL CONDITIONS #############################
// ###############################################################################
complex double *derivates = calloc(2*N_DGL,sizeof(double));
// Um den Anfangswert zu speichern, mit dem die Zwischen-Anfangswerte berechnet werden
complex double *temp = calloc(2*N_DGL,sizeof(double));
// Die errechneten Zwischenfunktionswerte - die Steigungen
complex double *k1 = calloc(2*N_DGL,sizeof(double));
complex double *k2 = calloc(2*N_DGL,sizeof(double));
        double *OUTPUT = calloc(100,sizeof(double));

complex double *deriv_stst = calloc(2*N_DGL,sizeof(double));

for (n=0;n<N_DGL;n++) derivates[n] = 0. + I* 0.;
for (n=0;n<100;n++)   OUTPUT[n]    =0.;
// set initial conditions: quantum in the ground state
rho(1,1) = 0. ; // QD in the Ground State with Zero Photons in the Cavity
rho(0,0) = 1. ; // QD in the Ground State with Zero Photons in the Cavity
rho(1,0) = 0. ; // QD in the Ground State with Zero Photons in the Cavity
rho(0,1) = 0. ; // QD in the Ground State with Zero Photons in the Cavity
OUTPUT[0] = 1.; // Trace=1
// ###########################################################################
// ########################## TIME _DOMAIN SOLUTION ##########################
// ###########################################################################
k=0; double progress_per_stst = 0.1; OUTPUT[99] = 100.;
// begin t integration

complex double W_emma = 0.;

double gs = 2.*(GAM/2.+gam_p);
double Gs = 2.*GAM;
double oms = 2.*OM;

double gp = (Gs+gs)/2.;
double gm = (Gs-gs)/2.;
double OMR = sqrt(oms*oms-Gs*Gs);

double W0 = rho(1,1) - rho(0,0);
double S0 = rho(1,1) + rho(0,0);
complex double B0 = rho(0,1) - rho(1,0);


while(OUTPUT[99]>0.0000000001)
{
 t=delta_t*k;

 W_emma  =  1.*(W0+(gs*Gs*S0)/((gp*gp)+(OMR*OMR)))*cexp(-gp*t)*cos(OMR*t);
 W_emma += -1.*((W0*gm) +(I*oms*B0) + (Gs*S0*(OMR*OMR+gm*gp))/((OMR*OMR)+(gp*gp)))*(cexp(-gp*t)*sin(OMR*t))/OMR;
 W_emma += -1.*(Gs*gs*S0)/((gp*gp)+(OMR*OMR));

 RUNGE_KUTTA(derivates,k1,k2,temp,OUTPUT,t,N_DGL);
 STEADY_STATE_TEST(derivates,deriv_stst,OUTPUT,Zahl_der_DGL);
 fprintf(f_spectrum,"%.10f\t %.10f\t %.10f %.10f \n",t,creal(rho(1,1)-rho(0,0)),2.*creal(W_emma),creal(rho(1,1)-rho(0,0)-W_emma));
 if (OUTPUT[99]<progress_per_stst)
    {
     progress_per_stst = progress_per_stst/3.;
     printf("%.0f - "      ,k*100./N);
     printf("rho(0,0)=%.5f - "   ,creal(rho(0,0)));
     printf("rho(1,1)=%.5f - "   ,creal(rho(1,1)));
     printf("ERROR=%.5f - ",creal(rho(1,1)-rho(0,0)-W_emma));
     printf("TRACE=%.10f - "    ,OUTPUT[0]);
     printf("StSt_TEST=%.10f \n",OUTPUT[99]);
    }
 k++;
} // end of time integration
// ########################## TAU_DOMAIN SOLUTION ############################
// free given arrays
free(temp);free(derivates);free(deriv_stst);free(k1);free(k2);free(OUTPUT);
fclose(f_spectrum);
return 0;
}

