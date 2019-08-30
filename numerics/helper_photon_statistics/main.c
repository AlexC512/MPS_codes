 #include <math.h>
 #include <stdlib.h>
 #include <complex.h>
 #include <stdio.h>

// ###################### Ende     ############################################  
///////////////////////////////////////////////////////////////////////////////  
// #########################  MAIN    ######################################### 
 int main() {  
int k;
double temp1,temp2,temp3,temp4,temp5,temp6,temp7, temp0,temp8,temp9,temp10,temp11,temp12;
temp1 = 0.; temp2 = 0.;

int N = 50;

FILE *f_dat;
f_dat = fopen("/home/alex/phot_stat_therm.dat","w");

double n_bar = 12.0;
double fak   = n_bar / (n_bar +1.);
temp2 = 1./(n_bar + 1.);

fprintf(f_dat," %.10f \t %.10f \n",0.,temp2);
temp1 = 1.;

for (k=1;k<N;k++)
{

temp1 = temp1*fak;
fprintf(f_dat," %.10f \t %.10f \n",1.*k,temp1*temp2);
}

fclose(f_dat);

FILE *f_dat2;
f_dat2 = fopen("/home/alex/phot_stat_coh.dat","w");

double alpha_quad = 12.0;
double fak_alpha   = exp(-alpha_quad);

fprintf(f_dat2," %.10f \t %.10f \n",0.,fak_alpha);
temp1 = 1.;

for (k=1;k<N;k++)
{

temp1 = temp1 * alpha_quad/(double)k;
fprintf(f_dat2," %.10f \t %.10f \n",1.*k,temp1*fak_alpha);
}

fclose(f_dat2);

FILE *f_dat3;
f_dat3 = fopen("/home/alex/phot_stat_fock.dat","w");

fprintf(f_dat3," %.10f \t %.10f \n",0.,0.);

for (k=1;k<N;k++)
{
temp1 = 0.;
fprintf(f_dat3," %.10f \t %.10f \n",1.*k,0.);
}

fclose(f_dat3);


  return 0;
}
