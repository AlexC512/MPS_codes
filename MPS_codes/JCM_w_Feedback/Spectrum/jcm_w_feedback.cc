#include "itensor/all.h"

using namespace itensor;

// #####################################################
void SETUP_MPS(MPS& psi,
               std::vector<Index>& bin,std::vector<Index>& binlink,
               int ttotal,int Nfb,int tls,int cav,int Ncav,
               Real init_EE,Real init_GG, Real init_pn[])
{
ITensor MPSTensor = ITensor();    
// the first Nfb psi(i) are loaded with empty bins for the first feedback loop
for(int j = 1; j<=Nfb; j++)
    {
      if ( j==1 ) 
      {    
      MPSTensor = ITensor(binlink[j],bin[j]);
      MPSTensor.set(binlink[j](1),bin[j](1),1.);
      }
      else
      {    
      MPSTensor = ITensor(binlink[j-1],bin[j],binlink[j]);
      MPSTensor.set(binlink[j-1](1),bin[j](1),binlink[j](1),1.);
      }
      psi.setA((j),MPSTensor);
    }
    
MPSTensor = ITensor(bin[tls],bin[cav],binlink[Nfb],binlink[Nfb+1]);
for(int phot=1;phot<Ncav;phot++)
    {    
    MPSTensor.set(bin[tls](2),bin[cav](phot),binlink[Nfb](1),binlink[Nfb+1](1),init_EE*init_pn[phot]);
    MPSTensor.set(bin[tls](1),bin[cav](phot),binlink[Nfb](1),binlink[Nfb+1](1),init_GG*init_pn[phot]);
    }
psi.setA(Nfb+1,MPSTensor);

for(int j = Nfb+2; j<=ttotal; j++)
    {
      MPSTensor = ITensor(binlink[j-1],bin[j-1],binlink[j]);
      MPSTensor.set(binlink[j-1](1),bin[j-1](1),binlink[j](1),1.);
      psi.setA((j),MPSTensor);
    }

MPSTensor = ITensor(binlink[ttotal],bin[ttotal]);
MPSTensor.set(binlink[ttotal](1),bin[ttotal](1),1.);
psi.setA((ttotal+1),MPSTensor);
// now I have an MPS, filled with bin[1] bin[2] ... bin[gtotal+1]
// set orthoCenter
psi.position(1);

return;
}
// #####################################################

// #####################################################
double emitter_population(ITensor A, int Nbin)
{
    Index s=findIndex(A,"sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(2),prime(s(2)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).real();
}
// #####################################################

// #####################################################
double cavity_population(ITensor A, int Ncav)
{
    Index s=findIndex(A,"cav") ;       
    ITensor cdgc = ITensor(s,prime(s));
    double cav_pop = 0.;
    
    for(int phot=1;phot<Ncav;phot++)
    cdgc.set(s(phot+1),prime(s(phot+1)),1.*phot);
    
    cav_pop = eltC( dag( A ) * noPrime(cdgc*A)).real();
    return cav_pop;
}
// #####################################################

// #####################################################
double bath_population(ITensor A, int Nbin)
{
    Index s=findIndex(A,"bath") ;       
    ITensor bdgb = ITensor(s,prime(s));
    for(int phot=1;phot<Nbin;phot++) 
    bdgb.set(s(phot+1),prime(s(phot+1)),1.*phot);
    return eltC(dag(prime(A,s))*bdgb*A).real();
}
// #####################################################


// #####################################################
Complex bath_coherence(ITensor A, int Nbin)
{
    Index s=findIndex(A,"bath") ;       
    ITensor SpSm = ITensor(s,prime(s));
    
    for(int phot=1;phot<Nbin;phot++) 
    SpSm.set(s(phot+1),prime(s(phot)),sqrt(phot*1.));
    
    return eltC( dag(A)*noPrime(SpSm*A) );
}
// #####################################################

// #####################################################
void SETUP_MPO(ITensor& U_evo,
               std::vector<Index>& bin,
               int Nfb, int tls,int cav, int Ncav,
               Complex coup_Ecav,Complex coup_gcav,Complex coup_bath)
{
double sqrt_phot;
auto H_sys  = ITensor(bin[cav],prime(bin[cav]));
for(int phot=1;phot<Ncav;phot++)
{    
 sqrt_phot=sqrt(1.*phot);   
 H_sys.set(bin[cav](phot  ),prime(bin[cav](phot+1)),coup_Ecav*sqrt_phot);
 H_sys.set(bin[cav](phot+1),prime(bin[cav](phot  )),coup_Ecav*sqrt_phot);
}
        
auto H_jcm  = ITensor(bin[tls],prime(bin[tls]),bin[cav],prime(bin[cav]));
for(int phot=1;phot<Ncav;phot++)
{    
 sqrt_phot=sqrt(1.*phot);   
 H_jcm.set(bin[tls](1),prime(bin[tls](2)),bin[cav](phot+1),prime(bin[cav](phot)  ),coup_gcav*sqrt_phot);
 H_jcm.set(bin[tls](2),prime(bin[tls](1)),bin[cav](phot)  ,prime(bin[cav](phot+1)),coup_gcav*sqrt_phot);    
} 
        
auto H_dis = ITensor(bin[cav],prime(bin[cav]),bin[Nfb+1],prime(bin[Nfb+1]));
// ATTENTION ... here I assume that the process of more photons in the bin as in the cavity is not 
// possible this is for sure correct, if no feedback is assumed
for(int phot=1;phot<Ncav;phot++)
{    
 for(int bphot=1;bphot<Ncav;bphot++)
 {    
  sqrt_phot=sqrt(1.*phot*bphot);      
  H_dis.set(bin[cav](phot+1),prime(bin[cav](phot  )),bin[Nfb+1](bphot  ),prime(bin[Nfb+1](bphot+1)), coup_bath*sqrt_phot);
  H_dis.set(bin[cav](phot  ),prime(bin[cav](phot+1)),bin[Nfb+1](bphot+1),prime(bin[Nfb+1](bphot  )),-coup_bath*sqrt_phot);
 }
}

auto H_fb = ITensor(bin[cav],prime(bin[cav]),bin[1],prime(bin[1]));
// ATTENTION ... here I assume that the process of more photons in the bin as in the cavity is not 
// possible this is for sure correct, if no feedback is assumed
for(int phot=1;phot<Ncav;phot++)
{    
 for(int bphot=1;bphot<Ncav;bphot++)
 {    
  sqrt_phot=sqrt(1.*phot*bphot);      
  H_fb.set(bin[cav](phot+1),prime(bin[cav](phot  )),bin[1](bphot  ),prime(bin[1](bphot+1)),-1.*coup_bath*sqrt_phot);
  H_fb.set(bin[cav](phot  ),prime(bin[cav](phot+1)),bin[1](bphot+1),prime(bin[1](bphot  )), 1.*coup_bath*sqrt_phot);
 }
}

// four indices fb-bin, bath-bin, tls bin, cavity-bin 
// Hsys[cav], H_jcm[tls,cav], H_dis[cav,bathbin],H_fb[cav,feedback]
H_sys = delta(bin[tls],prime(bin[tls])) * H_sys * delta(bin[1]    ,prime(bin[1])    ) ;
H_sys =                                   H_sys * delta(bin[Nfb+1],prime(bin[Nfb+1])) ;
H_jcm = delta(bin[1],prime(bin[1]))     * H_jcm * delta(bin[Nfb+1],prime(bin[Nfb+1])) ;
H_dis = delta(bin[tls],prime(bin[tls])) * H_dis * delta(bin[1],prime(bin[+1]))        ;
H_fb  = delta(bin[tls],prime(bin[tls])) * H_fb  * delta(bin[Nfb+1],prime(bin[Nfb+1])) ;      

auto H_int    =   H_sys 
                + H_jcm 
                + H_dis 
                + H_fb
              ;

auto H_int_1  =  H_int; 
auto H_int_2  = (1./2.)  * mapPrime(H_int_1*prime(H_int_1),2,1);
auto H_int_3  = (1./3.)  * mapPrime(H_int_1*prime(H_int_2),2,1);
auto H_int_4  = (1./4.)  * mapPrime(H_int_1*prime(H_int_3),2,1);
auto H_int_5  = (1./5.)  * mapPrime(H_int_1*prime(H_int_4),2,1);
auto H_int_6  = (1./6.)  * mapPrime(H_int_1*prime(H_int_5),2,1);
auto H_int_7  = (1./7.)  * mapPrime(H_int_1*prime(H_int_6),2,1);
auto H_int_8  = (1./8.)  * mapPrime(H_int_1*prime(H_int_7),2,1);
auto H_int_9  = (1./9.)  * mapPrime(H_int_1*prime(H_int_8),2,1);
auto H_int_10 = (1./10.) * mapPrime(H_int_1*prime(H_int_9),2,1);
 
auto temp_delta = ITensor(bin[tls],prime(bin[tls]));
temp_delta.set(bin[tls](1),prime(bin[tls](1)),1.);
temp_delta.set(bin[tls](2),prime(bin[tls](2)),1.);

temp_delta = temp_delta * delta(bin[cav],prime(bin[cav])); 
temp_delta = temp_delta * delta(bin[1],prime(bin[1])); 
temp_delta = temp_delta * delta(bin[Nfb+1],prime(bin[Nfb+1])); 

U_evo = temp_delta + H_int_1 + H_int_2 + H_int_3 + H_int_4 + H_int_5 
                   + H_int_6 + H_int_7 + H_int_8 + H_int_9 + H_int_10 
       ;
       
 return;      
}       
// ####################################################

// #####################################################
void SWAP_FORWARD(MPS& psi,const std::vector<Index>& bin, int from, int to, double cutoff)
{
    ITensor SWAP,U,S,V;
    Index iFB;
    for(int k=from;k<to;k++) // swap feedback bin next to tls bin
    {
            SWAP = psi.A(k)*psi.A(k+1); 
            iFB = findIndex(psi.A(k+1),"bath");
            if ( order(psi(k))==2 ) U=ITensor(iFB);
                              else  U=ITensor(iFB,commonIndex(psi.A(k),psi.A(k-1)));  
            svd(SWAP,U,S,V,{"Cutoff=",cutoff});
            psi.setA(k,U); 
            psi.setA(k+1,V*S); // orthoCenter wanders to system bin
    }
} 
// #####################################################

// #####################################################
void SWAP_BACKWARD(MPS& psi,const std::vector<Index>& bin, int from, int to, double cutoff)
{
    ITensor SWAP,U,S,V;
    Index iFB;
     for(int k=from;k>to;k--)
     {      
            SWAP = psi.A(k)*psi.A(k-1);
            iFB = findIndex(psi.A(k),"bath"); 
            if ( order(psi(k-1)) == 2 )  U=ITensor(iFB); // order(psi(k-1))==2, no commonIndex
                                  else   U=ITensor(iFB,commonIndex(psi.A(k-1),psi.A(k-2)));
            svd(SWAP,U,S,V,{"Cutoff=",cutoff}); 
            // -- lasse Orthocenter beim vorletzten Bin, der im naechsten Schritt wieder zum Systembin geswappt wird
            if (k-1==to)   { psi.setA(k-1,U);       psi.setA(k,V*S); }
            else           { psi.setA(k-1,U*S);     psi.setA(k,V);   }
     }
} 
// #####################################################

// #####################################################
void PRINT_MPS(MPS psi, int from, int to)
{
    for(int i=from;i<=to;i++) {printf("%i:",i); Print(psi(i)); }
    printf("------------------------------------ \n");
}
// #####################################################

// #####################################################
int main(int argc, char* argv[])
{    
// https://itensor.org/docs.cgi?page=tutorials/input    
if(argc != 2) { printfln("Usage: %s inputfile",argv[0]); return 0;     }
      
auto input = InputGroup(argv[1],"input");
//timestep
Real dt = input.getReal("time_step"); 
//end of integration
int t_end = input.getInt("time_end");
int Nfb   = input.getInt("feedback_at");
//coherent pumping strengths
Real E_cav = input.getReal("E_cav",0.); // if Omega_TLS is not given it return 0.
//TLS decay rate
Real Gamma = input.getReal("Gamma");
Real gcav  = input.getReal("gcav");
//dimension of local Hilbert space of each bin
int Ndim = input.getInt("Ndim");
int Nbin = input.getInt("Nbin"); 
int Ncav = input.getInt("Ncav"); 
//cutoff of schmidt values
Real cutoff = input.getReal("svdcutoff");
//maximal number of schmidtvalues 
int maxm = input.getInt("maxnumberofSV");
Real init_EE = input.getReal("init_EE");
Real init_GG = input.getReal("init_GG");
if ( fabs(init_EE+init_GG-1) > 0.00001 ) { printf("Emitter not properly initialized!!!\n"); return 0; }
Real init_pn[Ncav]; for(int pn=0;pn<=Ncav;pn++) init_pn[pn]=0.;
init_pn[1] = input.getReal("init_p0");
init_pn[2] = input.getReal("init_p1");
init_pn[3] = input.getReal("init_p2");
init_pn[4] = input.getReal("init_p3");
for(int pn=0;pn<=Ncav;pn++) init_pn[0] += init_pn[pn];
if (fabs(init_pn[0]-1.000)>0.00001) { printf("Cavity not initialized!!!\n"); return 0; }
int SHOW_EVERY_STEP = input.getInt("SHOW_EVERY_STEP");
int SHOW_EVERY_BATH_STEP = input.getInt("SHOW_EVERY_BATH_STEP");
int IFSPECTRUM = input.getInt("IFSPECTRUM");
int SpectrumSteps = input.getInt("SpectrumSteps");
Real SpectrumIntervall = input.getReal("SpectrumIntervall");
// --------------- setup output file -------------------
time_t curtime; struct tm *loctime; curtime = time(NULL); loctime = localtime (&curtime); 
int hour = loctime -> tm_hour;int minute = loctime -> tm_min;int second = loctime -> tm_sec;
int year = 1900+loctime -> tm_year;int month = loctime -> tm_mon;int day  = loctime -> tm_mday;
printf("Date: %d.%d.%.d -- Time: %d:%d:%d  \n",day,month+1,year,hour,minute,second);  

char FILE_NAME[2048+1024];
snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_System_Dynamics_w_Feedback_at_%i_Steps%i_Ncav_%i_dt_%.4f.dat",year,month+1,day,hour,minute,second,Nfb,t_end,Ncav,dt);
    
FILE *f_dat;
f_dat = fopen(FILE_NAME,"w");

// ###################################################################################    
// ###################################################################################    
// ###################################################################################
int ttotal = t_end+Nfb;
// setup physical degrees of freedom bins from 1 to gtotal+1
auto bin = std::vector<Index>(ttotal+5);
int tls = ttotal+3; 
int cav = ttotal+4;
bin.at(tls) = Index(Ndim,"sys");
bin.at(cav) = Index(Ncav,"cav");
for(int j = 1; j <= ttotal+1; ++j){ bin.at(j) = Index(Nbin,"bath"); }    

// setup links in between the physical degrees of freedom
auto binlink = std::vector<Index>(ttotal+2);
for(int j = 1; j <= ttotal+1; ++j){ binlink.at(j) = Index(1,"link");   }

MPS psi   = MPS(ttotal+1); 
SETUP_MPS(psi,bin,binlink,ttotal,Nfb,tls,cav,Ncav,init_EE,init_GG,init_pn);
//PRINT_MPS(psi,1,ttotal+1);
// ###################################################################################
// ###################################################################################
// ###################################################################################
Complex coup_Ecav = -Cplx_i*dt*E_cav;
Complex coup_gcav = -Cplx_i*dt*gcav;
Complex coup_bath =  sqrt(dt)*Gamma;

ITensor U_evo; 
SETUP_MPO(U_evo,bin,Nfb,tls,cav,Ncav,coup_Ecav,coup_gcav,coup_bath);    
//Print(U_evo);
// ###################################################################################
// ###################################################################################
// ###################################################################################
ITensor U,S,Vp,temp;
int j;
 
double rhoGE=0.;
double cavity_pop  = 1*init_pn[2]+2*init_pn[3]+3*init_pn[4];
double emitter_pop = init_EE;
double bath_pop    = 0.;
Complex bath_coh = 0.;
double sum_bath    = 0.;
double mps_norm = norm(psi(1)); 

Index iFB,iCB,iLink;

for(j = Nfb+1; j<=ttotal;j++)
{
 if( j%SHOW_EVERY_STEP == 0 ) 
   { 
    printf("Step: %i of  %i -- ",j-Nfb,t_end);
    printf("cav=%.10f -- em=%.10f -- bath=%.10f -- ",cavity_pop,emitter_pop,bath_pop);  
    printf("sum_bath=%.10f -- bath_coh=%.10f norm=%.10f",sum_bath,bath_coh,mps_norm );
    printf("\n");
   } 
 fprintf(f_dat,"%.10f \t",(j-Nfb)*dt); 
 fprintf(f_dat,"%.10f \t %.10f \t",cavity_pop,emitter_pop); 
 fprintf(f_dat,"%.10f \t %.10f \t",bath_coh.real(),bath_pop); 
 fprintf(f_dat,"\n"); fflush(f_dat);
 // ----------------------------------------------------------
 SWAP_FORWARD(psi,bin,j-Nfb,j-1,cutoff);
 // --------- apply mpo --------------------------------------
 iFB=findIndex(psi(j-1),"bath"); iCB=findIndex(psi(j+1),"bath");
 iLink=commonIndex(psi(j-2),psi(j-1));
 temp=noPrime(U_evo * psi(j-1)*psi(j)*psi(j+1)); 
 U=ITensor(iFB,iCB,iLink);
 svd(temp,U,S,Vp);

 psi.setA(j+1,Vp); 
 emitter_pop = emitter_population(Vp*S,Ndim);
 cavity_pop  = cavity_population(Vp*S,Ncav);
 mps_norm = eltC( dag(Vp*S)*Vp*S).real();
 
 temp=U*S;
 U=ITensor(iFB,iLink);
 svd(temp,U,S,Vp);
 
 psi.setA(j,Vp); 
 psi.setA(j-1,U*S); 
 bath_pop = bath_population(U*S,Nbin); sum_bath +=bath_pop;
 bath_coh = bath_coherence(U*S,Nbin);
 // ----------------------------------------------------------
 SWAP_BACKWARD(psi,bin,j-1,j-Nfb,cutoff);
 // ------------------prepare new mpo ------------------------
 if(j<ttotal)
 {    
 U_evo=delta(bin[j],bin[j+1])*U_evo*delta(prime(bin[j]),prime(bin[j+1])); 
 U_evo=delta(bin[j-Nfb],bin[j+1-Nfb])*U_evo*delta(prime(bin[j-Nfb]),prime(bin[j+1-Nfb]));
 }
 // Print(U_evo);  PRINT_MPS(psi,1,ttotal+1);
}

// ------------------- INFORMATION ABOUT CALCULATION LENGTH -------------------------------------------------
int ehour = loctime -> tm_hour; int eminute = loctime -> tm_min; int esecond = loctime -> tm_sec;
int eyear = 1900+loctime -> tm_year; int emonth = loctime -> tm_mon; int eday  = loctime -> tm_mday;
// give information how long it took
curtime = time(NULL); loctime = localtime (&curtime); 
ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
eyear = 1900+loctime -> tm_year; emonth = loctime -> tm_mon; eday  = loctime -> tm_mday;

printf("End of calculation: %d.%d.%.d -- Time: %d:%d:%d  \n",eday,emonth+1,eyear,ehour,eminute,esecond);  
printf("Time elapsed: %d.%d.%.d -- Time: %dh:%dmin:%dsec  -- in Minutes %d in Seconds %d  \n",eday-day,emonth-month,eyear-year,ehour-hour,eminute-minute,esecond-second,(ehour-hour)*60+(eminute-minute),(ehour-hour)*60*60+(eminute-minute)*60+(esecond-second));  

fprintf(f_dat,"# end of calculation %02d_%02d_%02d_%02d_%02d_%02d \n",eyear,emonth+1,eday,ehour,eminute,esecond);
fprintf(f_dat,"# time difference %02d_%02d_%02d_%02d_%02d_%02d -- in Minutes %d in Seconds %d \n",eyear-year,emonth-month,eday-day,ehour-hour,eminute-minute,esecond-second,(ehour-hour)*60+(eminute-minute),(ehour-hour)*60*60+(eminute-minute)*60+(esecond-second));
// ----------------- END INFORMATION ABOUT CALCULATION LENGTH -------------------------------------------------

printf("Steady-state = %.10f\n",4.*E_cav*E_cav/(Gamma*Gamma*Gamma*Gamma));
fclose(f_dat);

// ###################################################################################
// ###################################################################################
// ################################################################################### 

// ###################################################################################
// ###################################################################################
// ################################################################################### 

if (IFSPECTRUM == 1)
{    
    
snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_Bath_Dynamics_wo_Feedback_Steps%i_Ncav_%i_dt_%.4f.dat",year,month+1,day,hour,minute,second,t_end,Ncav,dt);
f_dat = fopen(FILE_NAME,"w");
// now we have a steady state and can calculate the spectrum <b^\dg(t)b(t-\tau)>
// setup the bath correlation array
auto bath_g1 = std::vector<Complex>(SpectrumSteps+1);
auto bath_G2= std::vector<Complex>(SpectrumSteps+1);
ITensor Tg1;
ITensor Tg2;
Index s,t;

// for feedback only take bins after one round trip into account
// ------------ swap ortoCenter to the first non-interacting bath bin -----------------
// in psi(ttotal+1) is the system // Print(psi(ttotal+1)); 
// in psi(ttotal-Nfb) is the first non-interacting bath bin ... // Print(iFB); Print(psi(ttotal-Nfb)); 
// in psi(ttotal-Nfb+1) is the orthoCenter // Print(psi(ttotal-Nfb+1));
j=ttotal-Nfb+1;
temp=psi(j)*psi(j-1); 
s=findIndex(psi(j-1),"bath"); //Print(s);
int RefBin=ttotal-Nfb;
U = ITensor(s,commonIndex(psi(j-2),psi(j-1))); //PrintData(U);
svd(temp,U,S,Tg1);
psi.setA(j-1,U*S);
psi.setA(j,Tg1); 
//Print(iFB);Print(psi(j-1));Print(psi(j)); 
// --------- end prepare via swap -------------------------------
// now define the flip operator for RefBin
ITensor Sp = ITensor(s,prime(s));
s=findIndex(psi(j-1),"bath");
for (int phot=1;phot<Nbin;phot++) Sp.set(s(phot),prime(s(phot+1)),sqrt(phot*1.)); //PrintData(Sp);
t=findIndex(psi(j-1),"bath"); 
ITensor Sm = ITensor(t,prime(t));
for (int phot=1;phot<Nbin;phot++) Sm.set(s(phot+1),prime(s(phot)),sqrt(phot*1.)); //PrintData(Sp);

// first we need the steady state value of the bath coherence
Complex bc_inf = bath_coherence(psi(j-1),Nbin);
double b_offset = abs(bc_inf)*abs(bc_inf);

Tg1= noPrime(Sm*psi(j-1)) ;
Tg1= dag( Tg1 ) * Tg1;
bath_g1[0]=eltC(Tg1);
double bdgb = eltC(Tg1).real(); 

Tg2= noPrime( Sm * noPrime(Sm*psi(j-1)));
bath_G2[0]=eltC( dag(Tg2)*Tg2 );

printf("b_c_inf=%.10f+I%.10f -- |b_c_inf|^2=%.10f -- <b^+b>=%.10f -- <b^+b^+bb>=%.10f\n",bc_inf.real(),bc_inf.imag(),b_offset,bdgb,bath_G2[0].real());

// normalize the g1 of the output-field
double bdiff = bath_g1[0].real()-b_offset;

for(int i=0;i<SpectrumSteps;i++)
{
// ---------------------------------------------------------------------------------------    
if ( i%SHOW_EVERY_BATH_STEP == 0) printf("%i: bath_g1=%.10f+I%.10f -- bath_G2=%.10f \n",i,(bath_g1[i].real()-b_offset)/bdiff,bath_g1[i].imag(),bath_G2[i].real()/(bdgb*bdgb));
fprintf(f_dat,"%.10f \t %.10f \t %.10f \t %.10f \t %.10f \n",i*dt,(bath_g1[i].real()-b_offset)/bdiff,bath_G2[i].real()/(bdgb*bdgb),bath_g1[i].real(),bath_g1[i].imag()); 
bath_g1[i]=(bath_g1[i].real()-b_offset)/bdiff; // for the spectrum the normalized coherence
// --------------------------------------------------------------------------------------

s=findIndex(psi(j-1-i),"bath"); 
Sp = ITensor(s,prime(s));
for (int phot=1;phot<Nbin;phot++) Sp.set(s(phot+1),prime(s(phot)),sqrt(1.*phot));

t=findIndex(psi(j-2-i),"bath"); 
Sm = ITensor(t,prime(t));
for (int phot=1;phot<Nbin;phot++) Sm.set(t(phot+1),prime(t(phot)),sqrt(1.*phot));
    
temp = psi(j-1-i)*psi(j-2-i); //Tg1=dag(temp)*temp;PrintData(Tg1);

Tg1= noPrime( Sp * dag(temp) ) * noPrime( Sm * temp);
bath_g1[i+1]=eltC(Tg1);

Tg2 = noPrime(Sp*noPrime(Sm*temp));
bath_G2[i+1]=eltC( dag(Tg2)*Tg2);

U = ITensor(s,commonIndex(psi(j-i-2),psi(j-i-3)));
svd(temp,U,S,Tg1);
psi.setA(j-i-2,U*S);
psi.setA(j-i-1,Tg1); 
}

fclose(f_dat);

// #######################################################################################################
// #######################################################################################################
// #######################################################################################################

snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_Spectrum_g2_wo_Feedback_Steps%i_Ncav_%i_dt_%.4f.dat",year,month+1,day,hour,minute,second,t_end,Ncav,dt);
f_dat = fopen(FILE_NAME,"w");

auto spectrum = std::vector<Complex>(SpectrumSteps+1);
double dw = SpectrumIntervall/(SpectrumSteps); 
double om;
   
for(int i=0;i<=SpectrumSteps-1;i++)
   {
     om = -0.5*(SpectrumIntervall)+i*dw;
     spectrum[i] = 0.;
     for(int b=0;b<=SpectrumSteps-1;b++)  spectrum[i] += exp(-Cplx_i*om*b*dt)*bath_g1[b];   
     fprintf(f_dat,"%.10f \t %.10f \t %.10f \n",om,spectrum[i].real(),spectrum[i].imag());
   }
fclose(f_dat);
// ifspectrum end */
} 
   
return 0;        
}
