#include "itensor/all.h"

using namespace itensor;


void MPO_SETUP(ITensor& U_evo, const std::vector<Index>& bin,int tls, int l_past, int l_now, int Nbin,Real Gamma_l, Real Gamma_r, Real phi_l, Real phi_r, Real dt        )
{
    
    int r_past= l_past+1;
    int r_now = l_now+1;
    
    // left atom interaction hamiltonian ... time local field (now)
    // define: |abc> == Atom left in a, Atom middle b, Atom right c, Atom 3, Atom 2, Atom 1, 
    Complex   gam_l    =  -Cplx_i*Gamma_l*sqrt(dt);
    Complex c_gam_l    =  -Cplx_i*Gamma_l*sqrt(dt);
    Complex   gam_l_fb =  -Cplx_i*Gamma_l*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_l);
    Complex c_gam_l_fb =  -Cplx_i*Gamma_l*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_l);

    Complex   gam_r    =  -Cplx_i*Gamma_r*sqrt(dt);
    Complex c_gam_r    =  -Cplx_i*Gamma_r*sqrt(dt);
    Complex   gam_r_fb =  -Cplx_i*Gamma_r*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_r);
    Complex c_gam_r_fb =  -Cplx_i*Gamma_r*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_r);
   
    auto H_dis_l  = ITensor(bin[tls],prime(bin[tls]),bin[l_now],prime(bin[l_now]));
    auto H_dis_r  = ITensor(bin[tls],prime(bin[tls]),bin[r_now],prime(bin[r_now]));
    double sqrt_phot=1.;
    for( int phot=1;phot<Nbin;phot++)
    {    
    sqrt_phot=sqrt(1.*phot);    
    // -------------------------------- left atom -----------------------------------------------------------------
    // for the right atom (2), local interactin with left moving photons -- we go from unprimed to primed    
    // <lr|, |lr> order in the initial states, 
    // in H_dis_l the right atom is loaded, because the right atom interacts first with left moving field
    H_dis_l.set(bin[tls](2),prime(bin[tls](1)),bin[l_now](phot),prime(bin[l_now](phot+1)),gam_r*sqrt_phot);
    H_dis_l.set(bin[tls](4),prime(bin[tls](3)),bin[l_now](phot),prime(bin[l_now](phot+1)),gam_r*sqrt_phot);        
    // for the left atom excitation means +4
    H_dis_l.set(bin[tls](1),prime(bin[tls](2)),bin[l_now](phot+1),prime(bin[l_now](phot)), c_gam_r*sqrt_phot);
    H_dis_l.set(bin[tls](3),prime(bin[tls](4)),bin[l_now](phot+1),prime(bin[l_now](phot)), c_gam_r*sqrt_phot);        
    // -------------------------------- right atom -----------------------------------------------------------------
    // for the left atom deexcitation means -4 we go from unprimed to primed    
    // in H_dis_r the left atom is loaded, because the left atom interacts first with the right moving field
    H_dis_r.set(bin[tls](3),prime(bin[tls](1)),bin[r_now](phot),prime(bin[r_now](phot+1)), gam_l*sqrt_phot);
    H_dis_r.set(bin[tls](4),prime(bin[tls](2)),bin[r_now](phot),prime(bin[r_now](phot+1)), gam_l*sqrt_phot);        
    // for the left atom excitation means +4
    H_dis_r.set(bin[tls](1),prime(bin[tls](3)),bin[r_now](phot+1),prime(bin[r_now](phot)), c_gam_l*sqrt_phot);
    H_dis_r.set(bin[tls](2),prime(bin[tls](4)),bin[r_now](phot+1),prime(bin[r_now](phot)), c_gam_l*sqrt_phot);        
    }
    
    // now bring dissipation hamiltonians into one product Hilbert space
    H_dis_l = delta(bin[r_now],prime(bin[r_now]))       * H_dis_l;    
    H_dis_r = delta(bin[l_now],prime(bin[l_now]))       * H_dis_r;
    
    auto H_dis =   H_dis_l + H_dis_r    ;  //PrintData(H_dis);
   
    // time non-local field with excitation exchange for left and right atom as the middle one has no extra feedback 
    auto H_fb_l  = ITensor(bin[tls],prime(bin[tls]),bin[r_past],prime(bin[r_past]));
    auto H_fb_r  = ITensor(bin[tls],prime(bin[tls]),bin[l_past],prime(bin[l_past]));
    
    for( int phot=1;phot<Nbin;phot++)
    {   
    sqrt_phot=sqrt(1.*phot);        
    H_fb_l.set(bin[tls](2),prime(bin[tls](1)),bin[r_past](phot),prime(bin[r_past](phot+1)), gam_r_fb*sqrt_phot);
    H_fb_l.set(bin[tls](4),prime(bin[tls](3)),bin[r_past](phot),prime(bin[r_past](phot+1)), gam_r_fb*sqrt_phot);        
    // the adjoint left atom interaction hamiltonian
    H_fb_l.set(bin[tls](1),prime(bin[tls](2)),bin[r_past](phot+1),prime(bin[r_past](phot)), c_gam_r_fb*sqrt_phot);
    H_fb_l.set(bin[tls](3),prime(bin[tls](4)),bin[r_past](phot+1),prime(bin[r_past](phot)), c_gam_r_fb*sqrt_phot);        
    // time non-local field ... 
    H_fb_r.set(bin[tls](3),prime(bin[tls](1)),bin[l_past](phot),prime(bin[l_past](phot+1)), gam_l_fb*sqrt_phot);
    H_fb_r.set(bin[tls](4),prime(bin[tls](2)),bin[l_past](phot),prime(bin[l_past](phot+1)), gam_l_fb*sqrt_phot);        
    // the adjoint left atom interaction hamiltonian
    H_fb_r.set(bin[tls](1),prime(bin[tls](3)),bin[l_past](phot+1),prime(bin[l_past](phot)), c_gam_l_fb*sqrt_phot);
    H_fb_r.set(bin[tls](2),prime(bin[tls](4)),bin[l_past](phot+1),prime(bin[l_past](phot)), c_gam_l_fb*sqrt_phot);        
    }

    auto H_fb  =   H_fb_l * delta(bin[l_past],prime(bin[l_past])) 
                 + H_fb_r * delta(bin[r_past],prime(bin[r_past]))                ;
                 
    H_fb  = delta(bin[l_now],prime(bin[l_now])) * H_fb * delta(bin[r_now],prime(bin[r_now]));
    H_dis = delta(bin[l_past],prime(bin[l_past])) * H_dis * delta(bin[r_past],prime(bin[r_past]));
           
    auto H_int    =  H_dis + H_fb;  //PrintData(H_int);           
    auto H_int_1  =  H_int; 
    
    auto temp_H_int_1 = prime(H_int_1);
    auto H_int_2  = (1./2.)  * mapPrime(temp_H_int_1*H_int_1,2,1); printf("2nd-order prepared: order of H=%i.\n",order(H_int_2));  
    auto H_int_3  = (1./3.)  * mapPrime(temp_H_int_1*H_int_2,2,1); printf("3rd-order prepared: order of H=%i.\n",order(H_int_3));
    auto H_int_4  = (1./4.)  * mapPrime(temp_H_int_1*H_int_3,2,1); printf("4th-order prepared: order of H=%i.\n",order(H_int_4));
    auto H_int_5  = (1./5.)  * mapPrime(temp_H_int_1*H_int_4,2,1); printf("5th-order prepared: order of H=%i.\n",order(H_int_5));
    auto H_int_6  = (1./6.)  * mapPrime(temp_H_int_1*H_int_5,2,1); printf("6th-order prepared: order of H=%i.\n",order(H_int_6));
    auto H_int_7  = (1./7.)  * mapPrime(temp_H_int_1*H_int_6,2,1); printf("7th-order prepared: order of H=%i.\n",order(H_int_7));
    auto H_int_8  = (1./8.)  * mapPrime(temp_H_int_1*H_int_7,2,1); printf("8th-order prepared: order of H=%i.\n",order(H_int_8));
    auto H_int_9  = (1./9.)  * mapPrime(temp_H_int_1*H_int_8,2,1); printf("9th-order prepared: order of H=%i.\n",order(H_int_9));
    auto H_int_10 = (1./10.) * mapPrime(temp_H_int_1*H_int_9,2,1); printf("10th-order prepared: order of H=%i.\n",order(H_int_10));

    ITensor delta_temp= ITensor(bin[tls],prime(bin[tls]));
    delta_temp.set(bin[tls](1),prime(bin[tls](1)),1.);
    delta_temp.set(bin[tls](2),prime(bin[tls](2)),1.);
    delta_temp.set(bin[tls](3),prime(bin[tls](3)),1.);
    delta_temp.set(bin[tls](4),prime(bin[tls](4)),1.);

    //ITensor temp;
    delta_temp = delta(bin[l_now],prime(bin[l_now]))*delta_temp*delta(bin[l_past],prime(bin[l_past]));
    delta_temp = delta(bin[r_now],prime(bin[r_now]))*delta_temp*delta(bin[r_past],prime(bin[r_past]));
    // now the taylor expansion of the U_evo

    U_evo =   delta_temp + H_int_1 + H_int_2 + H_int_3 
                        + H_int_4 + H_int_5 + H_int_6 + H_int_7  + H_int_8 + H_int_9 + H_int_10 ;
                   
    return;               
}


void MPS_SETUP(MPS& psi, const std::vector<Index>& bin, const std::vector<Index>& binlink, int ttotal, int Nfb, int tls, Real Init[])
{    
    // set up the mps 
    // dummy tensor we need
    ITensor MPSTensor = ITensor();
    // first MPS slot is filled with a bin which has no link to the left
     MPSTensor = ITensor(bin[1],binlink[1]);
     MPSTensor.set(bin[1](1),binlink[1](1),1.); // Level besetzt im TLS
     psi.setA((1),MPSTensor);

    for(int j = 2; j<=Nfb; j++)
    {
     MPSTensor = ITensor(bin[j],binlink[j-1],binlink[j]);
     MPSTensor.set(bin[j](1),binlink[j-1](1),binlink[j](1),1.); // Level besetzt im TLS
     psi.setA((j),MPSTensor);
    }   
    // 4*(0,1)+2*(0,1)+(0,1)+1
    MPSTensor = ITensor(bin[tls],binlink[Nfb],binlink[Nfb+1]);    
    MPSTensor.set(bin[tls](1),binlink[Nfb](1),binlink[Nfb+1](1),Init[1]); // 00 
    MPSTensor.set(bin[tls](2),binlink[Nfb](1),binlink[Nfb+1](1),Init[2]); // 01
    MPSTensor.set(bin[tls](3),binlink[Nfb](1),binlink[Nfb+1](1),Init[3]); // 10
    MPSTensor.set(bin[tls](4),binlink[Nfb](1),binlink[Nfb+1](1),Init[4]); // 11
    psi.setA((Nfb+1),MPSTensor); 

    for(int j = Nfb+1; j<=ttotal-1; j++)
    {
     MPSTensor = ITensor(bin[j],binlink[j],binlink[j+1]);
     MPSTensor.set(bin[j](1),binlink[j](1),binlink[j+1](1),1.); // Level besetzt im TLS
     psi.setA((j+1),MPSTensor);
    }
    // the last slot of MPS is filled with a tensor with only links to the left
    MPSTensor = ITensor(bin[ttotal],binlink[ttotal]);
    MPSTensor.set(bin[ttotal](1),binlink[ttotal](1),1.); // Level besetzt im TLS
    psi.setA((ttotal+1),MPSTensor);
    // now we have created an mps with (2*t_end time+Nfb) bins + 1 system bin located so
    // that feedback bins are initially empty to the left
    
//    for(int l=1;l<=ttotal+1;l++) {printf("%i:\n",l); PrintData(psi(l));}
}  

// -----------------------------------------------------------------------------------
// ------------------ INITIALIZATION PROCEDURES --------------------------------------
// -----------------------------------------------------------------------------------

double state_prop(ITensor A,int i)
{
    //search for Index in Tensor (should be system)
    Index s=findIndex(A, "sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(i),prime(s(i)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).real(); 
}


double re_coherence23(ITensor A)
{
    //search for Index in Tensor (should be system)
    Index s=findIndex(A, "sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(2),prime(s(3)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).real(); 
}


double TLS_norm(ITensor A)
{
    return eltC( dag(A)*A).real(); 
}

 
// moves right bin at "from" to "to" keeping OrthoCenter at to-1 in the left bin
void SWAP_RIGHT_FORWARD(MPS& psi, const std::vector<Index>& bin, int from, int to, double cutoff)
{
  
    ITensor SWAP,U,S,V;
    Index iLeft=findIndex(psi(from-1),"left");
    Index iRight=findIndex(psi(from),"right");
    
    for(int k=from;k<to;k++)
    {    
    SWAP = psi(k)*psi(k+1);
    U=ITensor(iRight,commonIndex(psi(k+1),psi(k+2)));
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});
    psi.setA(k,V*S); 
    psi.setA(k+1,U);
    
    SWAP = psi(k-1)*psi(k);
    U=ITensor(iLeft,commonIndex(psi(k),psi(k+1)));
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});
    psi.setA(k-1,V); 
    psi.setA(k,S*U);
    }

}    

// moves left bin at "from" to "to" keeping OrthoCenter at to-2 in the left bin
void SWAP_LEFT_BACKWARD(MPS& psi, const std::vector<Index>& bin, int from, int to, double cutoff)
{
  
    ITensor SWAP,U,S,V;
    Index iLeft=findIndex(psi(from),"left");
    Index iRight=findIndex(psi(from+1),"right");
    
    for(int k=from;k>to;k--)
    {    
    SWAP = psi(k)*psi(k-1); //PrintData(SWAP);
    if (order(psi(k-1))==2) U=ITensor(iLeft);
    else U=ITensor(iLeft,commonIndex(psi(k-1),psi(k-2)));
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});
    psi.setA(k,V*S); 
    psi.setA(k-1,U);
    
    SWAP = psi(k+1)*psi(k);
    U=ITensor(iRight,commonIndex(psi(k),psi(k-1)));
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});
    
    if (k-1==to) { psi.setA(k,U  ); psi.setA(k+1,V*S);}
    else         { psi.setA(k,U*S); psi.setA(k+1,V);  }
    }

}    

void PrintMPS(MPS psi,int from, int to)
{
    printf(" -------------------------------------- \n");
    for(int l=from;l<=to;l++) 
    {
        printf("%i:",l);
        Print(psi(l));
    }
    printf(" -------------------------------------- \n");
}    



int main(int argc, char* argv[])
{ 
if(argc != 2) { printfln("Usage: %s inputfile",argv[0]); return 0; } 
auto input = InputGroup(argv[1],"input");
// size of timestep
Real dt = input.getReal("time_step");
// maximum number of time steps
int t_end = input.getInt("time_end");
Real Gamma_l = input.getReal("Gamma_l",0.); //Gamma_l = Gamma_l*sqrt(dt);
Real Gamma_r = input.getReal("Gamma_r",0.); //Gamma_r = Gammm_r*sqrt(dt);
//dimension of local Hilbert space of each bin
int Nbin = input.getInt("Nbin",4); 
//cutoff of schmidt values
Real cutoff = input.getReal("svdcutoff");
//maximal number of schmidtvalues 
int maxm = input.getInt("maxnumberofSV");    
int fb  = input.getInt("feedback_time");
int Nfb = 2*fb; // size of array is twice because of left and right moving photons and one feedback loop 
int ttotal = 2*t_end+Nfb; printf("Total Number of time bins: %i\n",ttotal);
// ATTENTION at the moment I do not take phase in between the states into account
double Init[10];
Init[1] = input.getReal("init_00"); Init[2] = input.getReal("init_01"); 
Init[3] = input.getReal("init_10"); Init[4] = input.getReal("init_11"); 
    
Real phi_l = input.getReal("phi_l");      
Real phi_r = input.getReal("phi_r"); 
// ----------------------------------------------------------------------------------
// --------------------- PREPARE OUTPUT FILE ----------------------------------------
// ----------------------------------------------------------------------------------
time_t curtime; struct tm *loctime; curtime = time(NULL); loctime = localtime (&curtime); 
int hour = loctime -> tm_hour;int minute = loctime -> tm_min;int second = loctime -> tm_sec;
int year = 1900+loctime -> tm_year;int month = loctime -> tm_mon;int day  = loctime -> tm_mday;
printf("Date: %d.%d.%.d -- Time: %d:%d:%d  \n",day,month+1,year,hour,minute,second);  

double pop1=0.;
double pop2=0.;
double prop000 = 
// 00=1, 01=2 ,10=3,11=4
pop1 = Init[2]*Init[2]+Init[4]*Init[4]; // ATTENTION Init_xyz is real
pop2 = Init[3]*Init[3]+Init[4]*Init[4]; // ATTENTION Init_xyz is real

FILE *file;
FILE *f_prob;
char FILE_NAME[2048+1024];
snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_N=2_FB_at_Step=%i_of%i_phi_l_%.2f_phi_r_%.2f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,fb,t_end,phi_l,phi_r,dt,cutoff*1000000.);
file = fopen(FILE_NAME,"w");
fprintf(file,"## Date: %d.%d.%.d -- Time: %d:%d:%d  \n",day,month+1,year,hour,minute,second); 
fprintf(file,"## Gamma_l=%.10f - Gamma_m=%.10f - Gamma_r=%.10f \n",Gamma_l,0.,Gamma_r);    
fprintf(file,"## phi_l=%.10f - phi_m=%.10f - phi_r=%.10f \n",phi_l,0.,phi_r);    
fprintf(file,"## init_000=%.2f - init_001=%.2f - init_010=%.2f - init_011=%.2f \n",Init[1],Init[2],Init[3],Init[4]);    
fprintf(file,"## dt=%.6f - t_end=%i - fb=%i - cutoff=%.12f -- Nbin=%i \n",dt,t_end,fb,cutoff,Nbin);    

snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_N=2_Prob_FB_at_Step=%i_of%i_phi_l_%.2f_phi_r_%.2f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,fb,t_end,phi_l,phi_r,dt,cutoff*1000000.);
f_prob = fopen(FILE_NAME,"w");
fprintf(f_prob,"## Date: %d.%d.%.d -- Time: %d:%d:%d  \n",day,month+1,year,hour,minute,second); 
fprintf(f_prob,"## Gamma_l=%.10f - Gamma_m=%.10f - Gamma_r=%.10f \n",Gamma_l,0.,Gamma_r);    
fprintf(f_prob,"## phi_l=%.10f - phi_m=%.10f - phi_r=%.10f \n",phi_l,0.,phi_r);    
fprintf(f_prob,"## init_000=%.2f - init_001=%.2f - init_010=%.2f - init_011=%.2f \n",Init[1],Init[2],Init[3],Init[4]);    
fprintf(f_prob,"## dt=%.6f - t_end=%i - fb=%i - cutoff=%.12f -- Nbin=%i \n",dt,t_end,fb,cutoff,Nbin);    
// ----------------------------------------------------------------------------------
// --------------------- SETUP THE MPS ----------------------------------------------
// ----------------------------------------------------------------------------------
// create the physical indices
    
// array starts with zero ... so ttotal + 1 ... ttotal time bins, 1 system bin
auto bin = std::vector<Index>((int)ttotal+1); 
for(int j = 1; j <= ttotal; ++j){ if (j%2==1) bin.at(j) = Index(Nbin,"left");
                                      else        bin.at(j) = Index(Nbin,"right");         }
int tls=0; bin.at(tls) = Index(4,"sys"); 
// bin(0) system consists of three atoms 
auto binlink = std::vector<Index>((int)ttotal+1);
for(int j = 1; j <= ttotal; ++j) { binlink.at(j) = Index(1,"Link");  }
    
MPS psi=MPS(ttotal+1);
MPS_SETUP(psi,bin,binlink,ttotal,Nfb,tls,Init);

int sys_at=Nfb+1;
psi.position(1);
printf("Orthocenter at psi(%i) \n",orthoCenter(psi));

// ----------------------------------------------------------------------------------
// ---------------------END MPS SETUP -----------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ---------------------U_EVO SETUP -------------------------------------------------
// ----------------------------------------------------------------------------------  

// ATTENTION --- minimum feedback length between right and left emitter fb=4 bins 
// for fb=3 asymmetric feedback, and fb=2 no distance between left middle right
// for asymmetric feedback define left distance, right distance bin .... TODO

int l_past= sys_at-Nfb;
int r_past= sys_at-Nfb+1;
int l_now = Nfb+1; 
int r_now = Nfb+2;

ITensor U_evo; 
MPO_SETUP(U_evo,bin,tls,l_past,l_now,Nbin,Gamma_l,Gamma_r,phi_l,phi_r,dt);


// ----------------------------------------------------------------------------------
// --------------------- END U_EVO SETUP --------------------------------------------
// ----------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------
// --------------------- START EVOLUTION  -------------------------------------------
// ----------------------------------------------------------------------------------
double mps_norm=1.; mps_norm = TLS_norm(psi(sys_at));
Real output_r = 0.;
Real output_l = 0.;
ITensor BdgB,b_out;

Real Prob[10]; 
Prob[5]=0.;
for(int p=1;p<=4;p++) {Prob[p]=state_prop(psi(sys_at),p); Prob[5] += Prob[p];}
               
ITensor U,S,V,W,SWAP,MPO;
Index iFBl,iFBr,iCBl,iCBr; // index fuer feedback, system, und current bin 
Index iLinkp; 
for(int m=0;m<t_end;m++)
{   
    l_past= sys_at-Nfb    ;
    r_past= sys_at-Nfb+1  ;
    l_now = sys_at        ; 
    r_now = sys_at+1;
    // --- Status and File output ----------------------------------------------------------------------------------------
    printf("Step %i of %i: norm=%.10f -- ",m,t_end,mps_norm);
    printf("pop_l=%.10f -- pop_r=%.10f -- Sum_Prob=%.2f ",pop2,pop1,Prob[5]); 
    printf("output_l=%.10f -- output_r=%.10f",output_l,output_r); 
    printf("\n");
    fprintf(file,"%.10f  \t %.10f \t %.10f \t %.10f \t %.10f \t %.2f \n",m*dt,pop2,pop1,output_l,output_r,mps_norm);
    fflush(file); 
    fprintf(f_prob,"%.10f \t",m*dt); for(int p=1;p<=5;p++) fprintf(f_prob,"%.10f \t",Prob[p]); fprintf(f_prob,"\n");
    fflush(f_prob);
    // --- End - Status and File output -----------------------------------------------------------------------------------
    
    // -------------- PREPARE MPS TO APPLY MPO ------------------------------------------
    SWAP_RIGHT_FORWARD(psi,bin,r_past,sys_at-1,cutoff);   // move beyond middle bins to absorb orthoCenter
    // --------------------------------------------------------
    //PrintMPS(psi,1,ttotal+1);Print(U_evo);
    // --------------------------------------------------------
    // ----------- APPLY MPO START ----------------------------
    iFBl = findIndex(psi.A(sys_at-2),"left"); 
    iFBr = findIndex(psi.A(sys_at-1),"right");      
    iCBl = findIndex(psi.A(sys_at+1),"left");
    iCBr = findIndex(psi.A(sys_at+2),"right");    
    iLinkp = commonIndex(psi.A(sys_at-2),psi.A(sys_at-3));
    
    MPO = noPrime(psi(sys_at-2)*psi(sys_at-1)*U_evo*psi(sys_at)*psi(sys_at+1)*psi(sys_at+2));  
    U=ITensor(iLinkp,iFBl,iFBr,iCBl,iCBr); 
    svd(MPO,U,S,V,{"Cutoff=",cutoff});
    // tls in the next mps slot 
    psi.setA(sys_at+2,V); 
    Prob[5]=0;
    for(int p=1;p<=4;p++) {Prob[p]=state_prop(V*S,p); Prob[5] += Prob[p]; }

    pop1 =  Prob[2]+Prob[4];
    pop2 =  Prob[3]+Prob[4];
    mps_norm = eltC( dag(V*S)*V*S).real(); 
    // now factorize step by step the reservoir tensor moving orthoCenter
    MPO=U*S; U=ITensor(iLinkp,iFBl,iFBr,iCBl); svd(MPO,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at+1,V); // iCBr 
    MPO=U*S; U=ITensor(iLinkp,iFBl,iFBr); svd(MPO,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at,V); // iCBl     
    MPO=U*S; U=ITensor(iLinkp,iFBl); svd(MPO,U,S,V,{"Cutoff=",cutoff});  
    // -------- bath right -------------------------
    psi.setA(sys_at-1,V); // iFBr
    BdgB = ITensor(iFBr,prime(iFBr));
    BdgB.set(iFBr(2),prime(iFBr)(2),1.);
    b_out = V*S;
    output_r = eltC( dag( b_out )* noPrime(BdgB* b_out)  ).real(); 
    // --------------------------------------------
    // -------- bath right -------------------------
    psi.setA(sys_at-2,U*S); // iFBl 
    BdgB = ITensor(iFBl,prime(iFBl));
    BdgB.set(iFBl(2),prime(iFBl)(2),1.);
    b_out = U*S;
    output_l = eltC( dag( b_out )* noPrime(BdgB* b_out) ).real(); 
    // --------------------------------------------
    
    // ----------- SWAP BACK START -----------------------------
    SWAP_LEFT_BACKWARD(psi,bin,sys_at-2,l_past,cutoff);  
    // ----------- SWAP BACK END ------------------------------
    // --------------------------------------------------------
    // Update MPO, only if a further loop is possible
    if (m<t_end-1)
    {    
     U_evo = delta(bin[l_past],bin[l_past+2])*U_evo*delta(prime(bin[l_past]),prime(bin[l_past+2])) ;
     U_evo = delta(bin[r_past],bin[r_past+2])*U_evo*delta(prime(bin[r_past]),prime(bin[r_past+2])) ;
     U_evo = delta(bin[l_now],bin[l_now+2])*U_evo*delta(prime(bin[l_now]),prime(bin[l_now+2])) ;
     U_evo = delta(bin[r_now],bin[r_now+2])*U_evo*delta(prime(bin[r_now]),prime(bin[r_now+2])) ;
    }
    // increase the position of the system by 2!!
    sys_at++;sys_at++;
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

fprintf(file,"# end of calculation %02d_%02d_%02d_%02d_%02d_%02d \n",eyear,emonth+1,eday,ehour,eminute,esecond);
fprintf(file,"# time difference %02d_%02d_%02d_%02d_%02d_%02d -- in Minutes %d in Seconds %d \n",eyear-year,emonth-month,eday-day,ehour-hour,eminute-minute,esecond-second,(ehour-hour)*60+(eminute-minute),(ehour-hour)*60*60+(eminute-minute)*60+(esecond-second));

fprintf(f_prob,"# end of calculation %02d_%02d_%02d_%02d_%02d_%02d \n",eyear,emonth+1,eday,ehour,eminute,esecond);
fprintf(f_prob,"# time difference %02d_%02d_%02d_%02d_%02d_%02d -- in Minutes %d in Seconds %d \n",eyear-year,emonth-month,eday-day,ehour-hour,eminute-minute,esecond-second,(ehour-hour)*60+(eminute-minute),(ehour-hour)*60*60+(eminute-minute)*60+(esecond-second));
// ----------------- END INFORMATION ABOUT CALCULATION LENGTH -------------------------------------------------

fclose(file);fclose(f_prob);
 
return 0;
}
