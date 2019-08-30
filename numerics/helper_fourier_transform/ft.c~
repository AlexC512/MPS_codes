
 #include <math.h>
 #include <stdlib.h>
 #include <complex.h>
 #include <stdio.h>
 
 #define Pi 3.141592653589793238462643383279
 #define drittel 0.3333333333333333
 #define sechstel 0.1666666666666667
 #define N_DGL 5 
 #define Nt 75000
 #define N_ft 500
 #define delta_t 0.4
 #define M 0.001
 #define BAND_GAP_FREQUENCY 1.5
 #define GAMMA 0.01

 #define PULSE_FREQUENCY 1.5
 #define pulsfl 1.0 
 #define t_null 2000.   
 #define tau_parameter 500.  

 #define FT_INTERVALL 0.2

// ############ SWITCHES ###########################

//#define FREQUENCY_AFTER_TIME_SOLUTION
#define FREQUENCY_WHILE_TIME_SOLUTION 

//#define FT_OPT_AFTER_TIME

// #define SOURCE_SIN
// #define SOURCE_COS
// #define SOURCE_ROTDAMPING
// #define SOURCE_GAUSS
   #define SOURCE_RK

   #define ROTATING_FRAME

//  #define SVE_APPROX

// nur fuer SOURCE_RK

// ########### puls ########################################################
      double tnull=t_null;
      double tau=tau_parameter;
      double pulsflaeche=pulsfl;
    
 double puls (double t)
   {  double argument,temp, amplitude;
      tau=tau_parameter;
       
      amplitude=pulsflaeche*sqrt(Pi)/tau;
      argument=-(t-tnull)*(t-tnull)/(tau*tau);
      temp=exp(argument);
      temp *= amplitude;
      return temp;
       }

// ############### equations ###############################################
void calculate_next_time_step (complex double* derivates, complex double* derivates_out,  double t, double pulswert, double phase_pol, double phase_licht)
    {

    #ifdef ROTATING_FRAME
    derivates_out[0] =-GAMMA * derivates[0]
                      +I * M *  pulswert * phase_pol * phase_licht ;
    derivates_out[1] =0.;
    #endif

    return ;
    }

// ######################## Runge Kutta #######################################

void calculate_slope (complex double* derivates_initial, complex double* derivates_temp, complex double* derivates_in,
                        complex double* derivates_out, double Faktor, int SWITCH)
  {     int i;  
        for (i=0; i<N_DGL; i++)
		{// Anfangswert wird gesichert, fuer folgende Iterationen
                  if (SWITCH == 1) {derivates_temp[i]   = derivates_initial[i]; }
                 // Errechnete Steigung mit Schrittweite multipliziert 
                  derivates_in[i]         *= delta_t;
                 // Neuer Anfangswert wird berechnet, mit dem naechster Schritt gestartet wird
                  derivates_out[i]         = derivates_temp[i]   + derivates_in[i]; 
                 // Neuer Funktionswert wird in vier Schritten errechnet, gewichtet mit Faktor
		  derivates_initial[i]    += derivates_in[i]   * Faktor;
        	 }
       return; 
  } 


// ####################### Stoppuhr  #########################################
 #define   ESC             27
 #define   cls()           printf("%c[2J",ESC)
 #define   cursor(a,b)     printf("%c[%d;%dH",ESC,a,b)		 
		 
 struct timeval ts;		 
		 
 double stoppuhr ()
  { gettimeofday(&ts, NULL); 
    return (ts.tv_sec*1000)+(ts.tv_usec/1000); 
	  }


void restzeit_anzeigen(int k, double zeit_start)
{  double min,sek;
   double zeit_ende, v,zwischen_zeit,rest_zeit;
   int N = Nt;

//   cls();
//   cursor(14,13);
   zwischen_zeit = stoppuhr();
   v = k / (zwischen_zeit - zeit_start);
   rest_zeit = (N - k)/ v ;
   min = floor(rest_zeit*0.001/60.);
   sek = 0.001 * rest_zeit - 60. * floor(min);
   printf("%.2f % - Restzeit: %.0f min %.0f sek.\n ",k*100./N,min,sek);
 
 return;
}
// ################ Ende Stoppuhr ################################################
// ############## netcdf ###########################################################

#include <netcdf.h>
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define NDIMS 1
#define FILE_NAME "/home/alex/spektrum.nc"
// ###################### netcdf - Ende    ###################################  
// #########################  MAIN   ################################# 
 int main()
 {  int N,n,k;
    double temp1,temp2,temp3,temp4,temp5,temp6,temp7, temp0,temp8,temp9,temp10,temp11,temp12, pol_im_t, t_sc, faktor;

    double delta_k;
    complex double tm1,tm2,tm3,tm4,transformierte, transformierter_puls,tm7,tm9,tm6,tm5; 

#ifdef SOURCE_RK
    // ARRAYS fuer RK4 und fuer die Ergebnisse reservieren 
    // Ergebnis-Arrays
    complex double *derivates = calloc(2*N_DGL,sizeof(double));
    // Die errechneten Zwischenfunktionswerte - die Steigungen
    complex double *m = calloc(2*N_DGL,sizeof(double));
    // Die errechneten Zwischen-Anfangswerte
    complex double *step = calloc(2*N_DGL,sizeof(double));
    // Um den Anfangswert zu speichern, mit dem die Zwischen-Anfangswerte berechnet werden
    complex double *temp = calloc(2*N_DGL,sizeof(double));
#endif


    N = Nt;
    complex double fwert_t[N], fwert_p[N], fwert_puls[N]; 
    delta_k=2.*Pi/N;
    faktor = 1./N;
    temp1 = 0.; temp2 = 0.;

   double t0,t;
   #ifdef SVE_APPROX
   double delta=2.*Pi/(N-100);
   #else
   double delta=2.*Pi/N; 
   #endif
   double cos_freq = 1.5 * Pi;

   complex double ft_faktor;
   char c;
#ifdef SOURCE_RK
   complex double phase_pulse, phase_pol, phase_pulse_fak, phase_pol_fak;
   phase_pol = 1. + I *0.;
   phase_pulse = 1. + I *0.;
   phase_pol_fak   = cexp( I * BAND_GAP_FREQUENCY * 0.5 * delta_t);
   phase_pulse_fak = cexp(-I * PULSE_FREQUENCY    * 0.5 * delta_t);
#endif


#ifdef FT_OPT_AFTER_TIME

double ft_mitte_plus,ft_mitte_minus;
int ft_start,ft_ende;

ft_mitte_minus = (2.*Pi - ( BAND_GAP_FREQUENCY + 0.5 * FT_INTERVALL)  * delta_t)/delta;
ft_mitte_plus  = (2.*Pi - ( BAND_GAP_FREQUENCY - 0.5 * FT_INTERVALL)  * delta_t)/delta;

ft_start = (int) ft_mitte_minus;
ft_ende  = (int) ft_mitte_plus;
#endif


    double start,ende; start=stoppuhr();
/////////////////////////////////////////////////
// ############### NETCDF - GENERELL#############
    int retval;
    size_t storepos_write=0;
    size_t storepos_read=0;
//////////////////////////////////////////////////
// ################ Netcdf - Init - WRITE ########
        int ncid_write, x_dimid;  
        int varid_w0,varid_w1,varid_w2,varid_w3;
        int dimids[NDIMS],x;
        int NX=Nt;
      #ifdef FT_OPT_AFTER_TIME
           NX=ft_ende - ft_start;
      #endif
         
        // Zuerst Datei schreiben mit gegebener Anzahl von Eintraegen 
       if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &ncid_write))) ERR(retval);
       if ((retval = nc_def_dim(ncid_write, "x", NX, &x_dimid)))     ERR(retval);
       // Dimension errechnet und als integer-wert uebergeben
       dimids[0] = x_dimid;
       // Variablennamen definieren
      if ((retval = nc_def_var(ncid_write, "omega", NC_DOUBLE, NDIMS, dimids, &varid_w0))) ERR(retval);
      if ((retval = nc_def_var(ncid_write, "ft_reell", NC_DOUBLE, NDIMS, dimids, &varid_w1))) ERR(retval);
      if ((retval = nc_def_var(ncid_write, "ft_imag", NC_DOUBLE, NDIMS, dimids, &varid_w2))) ERR(retval);
      if ((retval = nc_enddef(ncid_write))) ERR(retval);
       // Zaehler der Schreibposition im NETCDF-File */
// ##################### NETCDF - Init - WRITE #############
// ######################## Programm #######################
   

#ifdef FREQUENCY_AFTER_TIME_SOLUTION
for (k=0; k<N; k++)
       {
       t = k * delta_t;
      #ifdef SOURCE_SIN 
        fwert_p[k]=sin(cos_freq * k * delta_t);   
      #endif
      #ifdef SOURCE_COS 
        fwert_p[k]=cos(cos_freq * k * delta_t);   
      #endif
      #ifdef SOURCE_ROTDAMPING 
       fwert_p[k]=cexp((I*cos_freq - 0.01) *k * delta_t);
      #endif
      #ifdef SOURCE_GAUSS
       fwert_p[k]=cos(cos_freq *k * delta_t) * puls(t);
      #endif 
      #ifdef SOURCE_RK
       #ifdef SVE_APPROX 
          fwert_t[k]= puls(t) ;
          fwert_p[k] = derivates[0] ;
       #else
          fwert_t[k]= puls(t)  * conj(phase_pulse); 
          //conj, weil das konjugierte in die BWGL eingeht, damit dreht das Feld aber nicht.
          fwert_p[k] = derivates[0]  * phase_pol;
       #endif
     
       // Step 1 -----
          calculate_next_time_step(derivates, m , t , puls(t), phase_pol, phase_pulse); 
          calculate_slope(derivates, temp, m, step, sechstel, 1);
       // Step 2 -----
          phase_pol *= phase_pol_fak;
          phase_pulse *= phase_pulse_fak;

          calculate_next_time_step(step, m, t+delta_t*0.5 , puls(t+delta_t*0.5), phase_pol, phase_pulse); 
          calculate_slope(derivates, temp, m, step, drittel, 2);
       // Step 3 -----
          calculate_next_time_step(step, m, t+delta_t*0.5 , puls(t+delta_t*0.5),phase_pol, phase_pulse); 
          calculate_slope(derivates, temp, m, step, drittel, 3);
       // Step 4 -----
          phase_pol *= phase_pol_fak;
          phase_pulse *= phase_pulse_fak; 

          calculate_next_time_step(step, m, t+delta_t , puls(t+delta_t),phase_pol, phase_pulse); 
          calculate_slope(derivates, temp, m, step, sechstel, 4);
       // --- neuer Funktionswert berechnet ---
       #endif
       }

tm1 = 1. + I* 0.;

#ifdef FT_OPT_AFTER_TIME

   for (n=ft_start; n<ft_ende; n++)
     {
       transformierte=0. + I * 0.;
       transformierter_puls =  0. + I * 0.;
       #ifdef SOURCE_RK
       temp1 = ( 2.*Pi - n * delta) * (1./delta_t); // weil die BANDGAP Frequency mit einhalb multipliziert wird
       #else
       temp1 = n * delta;
       #endif
       
       tm2 = 1. + 0. * I;
       tm1 = cexp(I*n*delta);

      for (k=0; k<N; k++)
       { 
        transformierte += fwert_p[k] * tm2 ;
        transformierter_puls += fwert_t[k] * tm2;
        tm2 = tm2 * tm1; 
        } 
          
       #ifdef SOURCE_RK
       temp4 = cabs(transformierter_puls);
       tm5 = transformierte * conj(transformierter_puls)  ;
 
       if (temp4 > 0.0000001) 
            { temp2 = creal(tm5)/temp4; 
              temp3= cimag(tm5)/temp4; }      
       else { 
              temp2 =0.; 
              temp3 = 0.;
            }

       #else      
       temp2 = faktor  * creal(transformierte);
       temp3 = faktor  * cimag(transformierte);       
       #endif
              
   
       if ((retval = nc_put_var1_double(ncid_write, varid_w0, &storepos_write, &temp1))) ERR(retval);
       if ((retval = nc_put_var1_double(ncid_write, varid_w1, &storepos_write, &temp2))) ERR(retval);
       if ((retval = nc_put_var1_double(ncid_write, varid_w2, &storepos_write, &temp3))) ERR(retval);

       // Zaehlerposition weiterzaehlen
       storepos_write++;
       // Faktor hochzaehlen
       // Restzeit - Spielerei
       if (n % 100 == 0) restzeit_anzeigen(n,start); 
       } 
#else
   for (n=0; n<N; n++)
     {
       transformierte=0. + I * 0.;
       transformierter_puls =  0. + I * 0.;
       #ifdef SOURCE_RK
       temp1 = ( 2.*Pi - n * delta) * (1./delta_t); // weil die BANDGAP Frequency mit einhalb multipliziert wird
       #else
       temp1 = n * delta;
       #endif
       tm2 = 1. + 0. * I;

       tm1 = cexp(I*n*delta);

       for (k=0; k<N; k++)
       { 
        transformierte += fwert_p[k] * tm2 ;
        transformierter_puls += fwert_t[k] * tm2;
        tm2 = tm2 * tm1; 
        } 
          
       #ifdef SOURCE_RK
       temp4 = cabs(transformierter_puls);
       tm5 = transformierte * conj(transformierter_puls)  ;
 
       if (temp4 > 0.0000001) 
            { temp2 = creal(tm5)/temp4; 
              temp3= cimag(tm5)/temp4; }      
       else { temp2 =0.; 
              temp3 = 0.; }

       #else      
       temp2 = faktor  * creal(transformierte);
       temp3 = faktor  * cimag(transformierte);       
       #endif
              
   
       if ((retval = nc_put_var1_double(ncid_write, varid_w0, &storepos_write, &temp1))) ERR(retval);
       if ((retval = nc_put_var1_double(ncid_write, varid_w1, &storepos_write, &temp2))) ERR(retval);
       if ((retval = nc_put_var1_double(ncid_write, varid_w2, &storepos_write, &temp3))) ERR(retval);

       // Zaehlerposition weiterzaehlen
       storepos_write++;
       // Faktor hochzaehlen
       // Restzeit - Spielerei
       if (n % 100 == 0) restzeit_anzeigen(n,start); 
       } 
#endif

#endif

// nun die Zeitlösung ungespeichert, direkt in die Frequenzraumlösung transformieren

#ifdef FREQUENCY_WHILE_TIME_SOLUTION

complex double *ft_pol  = calloc(2*N,sizeof(double));
complex double *ft_puls = calloc(2*N,sizeof(double));
double pulswert;

for (n=0; n<N; n++)
       {
       t = n * delta_t;
       tm2 = 1. + 0. * I;
       tm1 = cexp(I*n*delta);
       pulswert = puls(t);
       #ifdef SOURCE_RK
       #ifdef SVE_APPROX 
       for (k=0; k<N; k++)
       { 
        ft_pol[k]  += derivates[0] * tm2 ;
        ft_puls[k] += pulswert * tm2;
        tm2 = tm2 * tm1; 
        } 
       #else
       for (k=0; k<N; k++)
        { 
          ft_pol[k]  += derivates[0]  * phase_pol * tm2 ;
          ft_puls[k] += pulswert  * conj(phase_pulse) * tm2;
          tm2 = tm2 * tm1; 
         } 
       #endif
     
       // Step 1 -----
          calculate_next_time_step(derivates, m , t , puls(t), phase_pol, phase_pulse); 
          calculate_slope(derivates, temp, m, step, sechstel, 1);
       // Step 2 -----
          phase_pol *= phase_pol_fak;
          phase_pulse *= phase_pulse_fak;

          calculate_next_time_step(step, m, t+delta_t*0.5 , puls(t+delta_t*0.5), phase_pol, phase_pulse); 
          calculate_slope(derivates, temp, m, step, drittel, 2);
       // Step 3 -----
          calculate_next_time_step(step, m, t+delta_t*0.5 , puls(t+delta_t*0.5),phase_pol, phase_pulse); 
          calculate_slope(derivates, temp, m, step, drittel, 3);
       // Step 4 -----
          phase_pol *= phase_pol_fak;
          phase_pulse *= phase_pulse_fak; 

          calculate_next_time_step(step, m, t+delta_t , puls(t+delta_t),phase_pol, phase_pulse); 
          calculate_slope(derivates, temp, m, step, sechstel, 4);
       // --- neuer Funktionswert berechnet ---
       #endif
       // Restzeit - Spielerei
       if (n % 100 == 0) restzeit_anzeigen(n,start); 

} 
 

  for (n=0; n<N; n++)
   {  
       #ifdef SOURCE_RK
       temp1 = ( 2.*Pi - n * delta) * (1./delta_t); // weil die BANDGAP_FREQ mit 0.5 multipliziert wird
       #else
       temp1 = n * delta;
       #endif

       #ifdef SOURCE_RK
       tm4   = ft_puls[n];  
       temp4 = cabs(tm1);
       tm5   = ft_pol[n] * conj(tm4)  ;
 
       if (temp4 > 0.0000001) 
            { temp2 = creal(tm5)/temp4; 
              temp3 = cimag(tm5)/temp4; }      
       else { temp2 =0.; 
              temp3 = 0.; }

       #else      
       temp2 = 0.2;
       temp3 = 0.2;       
       #endif
              
   
       if ((retval = nc_put_var1_double(ncid_write, varid_w0, &storepos_write, &temp1))) ERR(retval);
       if ((retval = nc_put_var1_double(ncid_write, varid_w1, &storepos_write, &temp2))) ERR(retval);
       if ((retval = nc_put_var1_double(ncid_write, varid_w2, &storepos_write, &temp3))) ERR(retval);

       // Zaehlerposition weiterzaehlen
       storepos_write++;
       // Faktor hochzaehlen
       } 

free(ft_pol);
free(ft_puls);
#endif
    ende=stoppuhr();
    printf("\n Rechenzeit: %.3f s\n",(ende-start)*0.001);
//    printf("\n FT_Mitte %i \n", ft_ende - ft_start  );
#ifndef SOURCE_RK
    printf("%.5f Resonanzfrequenz",delta_t * cos_freq);
#endif
  free(m);
  free(step);
  free(temp);
  free(derivates);
#error FALSCHER RUNGE KUTTA!!!   

 // ######### Netcdf - Datei schließen ##################
  
  if ((retval = nc_close(ncid_write)))  ERR(retval);
  // ##################################################### 
  return 0;
}
