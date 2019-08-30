 #include <math.h>
 #include <stdlib.h>
 #include <complex.h>
 #include <stdio.h>


#define INPUT_FILE "bath_corr.dat" //LENGH
#define OUTPUT_NAME "bath_corr_smooth.dat"

#define SpectrumSteps 25000
#define SpectrumIntervall 3.0

// ###################### Ende     ############################################  
///////////////////////////////////////////////////////////////////////////////  
// #########################  MAIN    ######################################### 
 int main() {  
int k;
double temp1,temp2,temp3,temp4,temp5,temp6,temp7, temp0,temp8,temp9,temp10,temp11,temp12;
temp1 = 0.; temp2 = 0.;

double abcissa[SpectrumSteps];
double ordinata[SpectrumSteps];  
temp1 = 0.;
temp2 = 0.;
temp3 = 0.;
temp4 = 0.;
temp5 = 0.;

FILE *f_source;
f_source=fopen(INPUT_FILE,"r");

k=0;

// find average offset

double offset=0.;

while (feof(f_source)==0) 
{
 fscanf(f_source,"%lf \t %lf \t %lf \t %lf \t %lf \n",&temp1,&temp2,&temp3,&temp4,&temp5);
 abcissa[k] = temp1;
 ordinata[k] = temp2;
 if ( (k>SpectrumSteps-1000) && (k<SpectrumSteps-200) ) offset += temp2/800.;
 k++;
}
fclose(f_source);

FILE *f_dat;
f_dat = fopen(OUTPUT_NAME,"w");
double RAD_DEPH = 0.0001250;

for (k=0; k<SpectrumSteps; k++) 
 {
  temp1=abcissa[k];
  ordinata[k] = (ordinata[k]-offset)*exp(-RAD_DEPH*temp1);
  fprintf(f_dat,"%.10f \t %.10f \n",temp1,ordinata[k] );
 }
fclose(f_dat);


complex double spectrum[SpectrumSteps];
double dw = SpectrumIntervall/(SpectrumSteps); 
double om;
double dt = abcissa[99]-abcissa[98];


f_dat = fopen("Spectrum.dat","w");

for(int i=0;i<=SpectrumSteps-1;i++)
   {
     om = -0.5*(SpectrumIntervall)+i*dw;
     spectrum[i] = 0.;
     for(int b=0;b<=SpectrumSteps-1;b++) 
        { 
                spectrum[i] += cexp(-I*om*b*dt)* ordinata[b];
        }
    fprintf(f_dat,"%.10f \t %.10f \t %.10f \n",om,creal(spectrum[i]),cimag(spectrum[i]));
   }
fclose(f_dat);


  return 0;
}
