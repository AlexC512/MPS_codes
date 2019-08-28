#include "itensor/all.h"

using namespace itensor;

double state4_pop(ITensor A)
{
    //search for Index in Tensor (should be system)
    Index s=findIndex(A, "sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(4),prime(s(4)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).real(); 
}

double state3_pop(ITensor A)
{
    //search for Index in Tensor (should be system)
    Index s=findIndex(A, "sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(3),prime(s(3)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).real(); 
}

double state2_pop(ITensor A)
{
    //search for Index in Tensor (should be system)
    Index s=findIndex(A, "sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(2),prime(s(2)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).real(); 
}

double state1_pop(ITensor A)
{
    //search for Index in Tensor (should be system)
    Index s=findIndex(A, "sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(1),prime(s(1)),1.);
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

    
    
int main(int argc, char* argv[])
    { 
               
    if(argc != 2) 
      { 
      //reminds us to give an input file if we forget
      printfln("Usage: %s inputfile",argv[0]); 
      return 0; 
      }
      
    auto input = InputGroup(argv[1],"input");
    // size of timestep
    Real dt = input.getReal("time_step");
    // maximum number of time steps
    int t_end = input.getInt("time_end");
    //TLS 1 decay rate
    Real Gamma1 = input.getReal("Gamma1",0.);
    //TLS 2 decay rate
    Real Gamma2 = input.getReal("Gamma2",0.);
    //dimension of local Hilbert space of each bin
    int Nbin = input.getInt("Nbin",4); 
    //cutoff of schmidt values
    Real cutoff = input.getReal("svdcutoff");
    //maximal number of schmidtvalues 
    int maxm = input.getInt("maxnumberofSV");    
    int fb  = input.getInt("feedback_time");
    int Nfb = 2*fb;
    // size of array is twice because of left and right moving photons and one feedback loop 
    int ttotal = 2*t_end+Nfb; printf("Total Number of time bins: %i\n",ttotal);
    
    Real init_EE   = input.getReal("init_EE"); 
    Real init_SYM  = input.getReal("init_SYM");
    Real init_ASYM = input.getReal("init_ASYM");
    Real init_GG   = input.getReal("init_GG");
    Real phi = input.getReal("phi"); 
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // --------------------- SETUP THE MPS ----------------------------------------------
    // ----------------------------------------------------------------------------------
    // create the physical indices
    
    // array starts with zero ... so ttotal + 1 ... ttotal time bins, 1 system bin
    auto bin = std::vector<Index>((int)ttotal+1); 
    for(int j = 1; j <= ttotal; ++j){ if (j%2==1) bin.at(j) = Index(3,"left");
                                        else bin.at(j) = Index(3,"right");         }
    int tls=0;
    bin.at(tls) = Index(4,"sys"); // bin(0) system consists of two atoms 
    // increase number by one because we do not use zeroth link!!
    auto binlink = std::vector<Index>((int)ttotal+1);
    for(int j = 1; j <= ttotal; ++j) { binlink.at(j) = Index(1,"Link");  }
    // we need space in the mps ttotal time bins + 1 system bin ... actually there is psi(0)
    // but it does not count?
    MPS psi=MPS(ttotal+1);
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
    // at MPS slot Nfb the tls tensor is put
    // I decide here whether I stay in the single excitation regime ATTENTION
    // otherwise the evolution operator must be increased in size!!! ATTENTION
    MPSTensor = ITensor(bin[tls],binlink[Nfb],binlink[Nfb+1]);
    MPSTensor.set(bin[tls](1),binlink[Nfb](1),binlink[Nfb+1](1),init_GG); 
    MPSTensor.set(bin[tls](2),binlink[Nfb](1),binlink[Nfb+1](1),init_SYM);
    MPSTensor.set(bin[tls](3),binlink[Nfb](1),binlink[Nfb+1](1),init_ASYM);
    MPSTensor.set(bin[tls](4),binlink[Nfb](1),binlink[Nfb+1](1),init_EE);
    int sys_at=Nfb+1;
    psi.setA((sys_at),MPSTensor); 

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
    
    psi.position(sys_at);
    printf("Orthocenter at psi(%i) \n",orthoCenter(psi));

//    int l; for(l=1;l<=ttotal+1;l++) {printf("%i:\n",l); PrintData(psi(l));}
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ---------------------END MPS SETUP -----------------------------------------------
    // ----------------------------------------------------------------------------------
  
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ---------------------U_EVO SETUP -------------------------------------------------
    // ----------------------------------------------------------------------------------
    int l_past= sys_at-Nfb;
    int r_past= sys_at-Nfb+1;
    int l_now = Nfb+1; 
    int r_now = Nfb+2;
    
    // stay in the single excitaton regime ... ATTENTION
    // left atom interaction hamiltonian ... time local field (now)
    Complex   gam_l    =  Gamma1*sqrt(dt);
    Complex c_gam_l    = -Gamma1*sqrt(dt);
    Complex   gam_l_fb = -gam_l*exp(-Cplx_i*phi);
    Complex c_gam_l_fb =  gam_l*exp( Cplx_i*phi);

    Complex   gam_r    =  Gamma2*sqrt(dt);
    Complex c_gam_r    = -Gamma2*sqrt(dt);
    Complex   gam_r_fb = -gam_r*exp(-Cplx_i*phi);
    Complex c_gam_r_fb =  gam_r*exp( Cplx_i*phi);
   
    auto H_dis_l = ITensor(bin[tls],prime(bin[tls]),bin[l_now],prime(bin[l_now]));
    auto H_dis_r = ITensor(bin[tls],prime(bin[tls]),bin[r_now],prime(bin[r_now]));
    for( int phot=1;phot<3;phot++)
    {    
    H_dis_l.set(bin[tls](2),prime(bin[tls](1)),bin[l_now](phot),prime(bin[l_now](phot+1)), gam_l*sqrt(1.*phot));
    H_dis_l.set(bin[tls](3),prime(bin[tls](1)),bin[l_now](phot),prime(bin[l_now](phot+1)),-gam_l*sqrt(1.*phot));        
    H_dis_l.set(bin[tls](4),prime(bin[tls](2)),bin[l_now](phot),prime(bin[l_now](phot+1)), gam_l*sqrt(1.*phot));        
    H_dis_l.set(bin[tls](4),prime(bin[tls](3)),bin[l_now](phot),prime(bin[l_now](phot+1)), gam_l*sqrt(1.*phot));        
    // the adjoint left atom interaction hamiltonian
    H_dis_l.set(bin[tls](1),prime(bin[tls](2)),bin[l_now](phot+1),prime(bin[l_now](phot)), c_gam_l*sqrt(1.*phot));
    H_dis_l.set(bin[tls](1),prime(bin[tls](3)),bin[l_now](phot+1),prime(bin[l_now](phot)),-c_gam_l*sqrt(1.*phot));        
    H_dis_l.set(bin[tls](2),prime(bin[tls](4)),bin[l_now](phot+1),prime(bin[l_now](phot)), c_gam_l*sqrt(1.*phot));        
    H_dis_l.set(bin[tls](3),prime(bin[tls](4)),bin[l_now](phot+1),prime(bin[l_now](phot)), c_gam_l*sqrt(1.*phot));        
    // right atom interation hamiltonian ... photons in time bins are created
    H_dis_r.set(bin[tls](2),prime(bin[tls](1)),bin[r_now](phot),prime(bin[r_now](phot+1)), gam_r*sqrt(1.*phot));
    H_dis_r.set(bin[tls](3),prime(bin[tls](1)),bin[r_now](phot),prime(bin[r_now](phot+1)), gam_r*sqrt(1.*phot));        
    H_dis_r.set(bin[tls](4),prime(bin[tls](2)),bin[r_now](phot),prime(bin[r_now](phot+1)), gam_r*sqrt(1.*phot));        
    H_dis_r.set(bin[tls](4),prime(bin[tls](3)),bin[r_now](phot),prime(bin[r_now](phot+1)),-gam_r*sqrt(1.*phot));        
    // the adjoint right atom interaction hamiltonian ... photons in time bins are annihilated
    H_dis_r.set(bin[tls](1),prime(bin[tls](2)),bin[r_now](phot+1),prime(bin[r_now](phot)), c_gam_r*sqrt(1.*phot));
    H_dis_r.set(bin[tls](1),prime(bin[tls](3)),bin[r_now](phot+1),prime(bin[r_now](phot)), c_gam_r*sqrt(1.*phot));        
    H_dis_r.set(bin[tls](2),prime(bin[tls](4)),bin[r_now](phot+1),prime(bin[r_now](phot)), c_gam_r*sqrt(1.*phot));        
    H_dis_r.set(bin[tls](3),prime(bin[tls](4)),bin[r_now](phot+1),prime(bin[r_now](phot)),-c_gam_r*sqrt(1.*phot));        
    }
    auto H_dis =   H_dis_l * delta(bin[r_now],prime(bin[r_now])) 
                 + H_dis_r * delta(bin[l_now],prime(bin[l_now])) 
               ;  


    /* // if every emitter has its own local and nonlocal field
    // interaction with the pastbin
    auto H_fb_l  = ITensor(bin[tls],prime(bin[tls]),bin[l_past],prime(bin[l_past]));
    auto H_fb_r  = ITensor(bin[tls],prime(bin[tls]),bin[r_past],prime(bin[r_past]));

    // time non-local field ... ATTENTION ... single excitation limit 
    H_fb_l.set(bin[tls](2),prime(bin[tls](1)),bin[l_past](1),prime(bin[l_past](2)), gam_l_fb);
    H_fb_l.set(bin[tls](3),prime(bin[tls](1)),bin[l_past](1),prime(bin[l_past](2)),-gam_l_fb);        
    H_fb_l.set(bin[tls](4),prime(bin[tls](2)),bin[l_past](1),prime(bin[l_past](2)), gam_l_fb);        
    H_fb_l.set(bin[tls](4),prime(bin[tls](3)),bin[l_past](1),prime(bin[l_past](2)), gam_l_fb);        
    // the adjoint left atom interaction hamiltonian
    H_fb_l.set(bin[tls](1),prime(bin[tls](2)),bin[l_past](2),prime(bin[l_past](1)), c_gam_l_fb);
    H_fb_l.set(bin[tls](1),prime(bin[tls](3)),bin[l_past](2),prime(bin[l_past](1)),-c_gam_l_fb);        
    H_fb_l.set(bin[tls](2),prime(bin[tls](4)),bin[l_past](2),prime(bin[l_past](1)), c_gam_l_fb);        
    H_fb_l.set(bin[tls](3),prime(bin[tls](4)),bin[l_past](2),prime(bin[l_past](1)), c_gam_l_fb);                
    // time non-local field ... ATTENTION ... single excitation limit 
    H_fb_r.set(bin[tls](2),prime(bin[tls](1)),bin[r_past](1),prime(bin[r_past](2)), gam_r_fb);
    H_fb_r.set(bin[tls](3),prime(bin[tls](1)),bin[r_past](1),prime(bin[r_past](2)), gam_r_fb);        
    H_fb_r.set(bin[tls](4),prime(bin[tls](2)),bin[r_past](1),prime(bin[r_past](2)), gam_r_fb);        
    H_fb_r.set(bin[tls](4),prime(bin[tls](3)),bin[r_past](1),prime(bin[r_past](2)),-gam_r_fb);        
    // the adjoint left atom interaction hamiltonian
    H_fb_r.set(bin[tls](1),prime(bin[tls](2)),bin[r_past](2),prime(bin[r_past](1)), c_gam_r_fb);
    H_fb_r.set(bin[tls](1),prime(bin[tls](3)),bin[r_past](2),prime(bin[r_past](1)), c_gam_r_fb);        
    H_fb_r.set(bin[tls](2),prime(bin[tls](4)),bin[r_past](2),prime(bin[r_past](1)), c_gam_r_fb);        
    H_fb_r.set(bin[tls](3),prime(bin[tls](4)),bin[r_past](2),prime(bin[r_past](1)),-c_gam_r_fb);        

    auto H_fb  =   H_fb_l * delta(bin[r_past],prime(bin[r_past])) 
                 + H_fb_r * delta(bin[l_past],prime(bin[l_past]))
               ;  
    */

    // time non-local field with excitation exchange ... ATTENTION ... single excitation limit 
    // interaction with the pastbin
    auto H_fb_l  = ITensor(bin[tls],prime(bin[tls]),bin[r_past],prime(bin[r_past]));
    auto H_fb_r  = ITensor(bin[tls],prime(bin[tls]),bin[l_past],prime(bin[l_past]));
    for( int phot=1;phot<3;phot++)
    {    
    H_fb_l.set(bin[tls](2),prime(bin[tls](1)),bin[r_past](phot),prime(bin[r_past](phot+1)), gam_l_fb*sqrt(1.*phot));
    H_fb_l.set(bin[tls](3),prime(bin[tls](1)),bin[r_past](phot),prime(bin[r_past](phot+1)),-gam_l_fb*sqrt(1.*phot));        
    H_fb_l.set(bin[tls](4),prime(bin[tls](2)),bin[r_past](phot),prime(bin[r_past](phot+1)), gam_l_fb*sqrt(1.*phot));        
    H_fb_l.set(bin[tls](4),prime(bin[tls](3)),bin[r_past](phot),prime(bin[r_past](phot+1)), gam_l_fb);        
    // the adjoint left atom interaction hamiltonian
    H_fb_l.set(bin[tls](1),prime(bin[tls](2)),bin[r_past](phot+1),prime(bin[r_past](phot)), c_gam_l_fb*sqrt(1.*phot));
    H_fb_l.set(bin[tls](1),prime(bin[tls](3)),bin[r_past](phot+1),prime(bin[r_past](phot)),-c_gam_l_fb*sqrt(1.*phot));        
    H_fb_l.set(bin[tls](2),prime(bin[tls](4)),bin[r_past](phot+1),prime(bin[r_past](phot)), c_gam_l_fb*sqrt(1.*phot));        
    H_fb_l.set(bin[tls](3),prime(bin[tls](4)),bin[r_past](phot+1),prime(bin[r_past](phot)), c_gam_l_fb*sqrt(1.*phot));                
    // time non-local field ... ATTENTION ... single excitation limit 
    H_fb_r.set(bin[tls](2),prime(bin[tls](1)),bin[l_past](phot),prime(bin[l_past](phot+1)), gam_r_fb*sqrt(1.*phot));
    H_fb_r.set(bin[tls](3),prime(bin[tls](1)),bin[l_past](phot),prime(bin[l_past](phot+1)), gam_r_fb*sqrt(1.*phot));        
    H_fb_r.set(bin[tls](4),prime(bin[tls](2)),bin[l_past](phot),prime(bin[l_past](phot+1)), gam_r_fb*sqrt(1.*phot));        
    H_fb_r.set(bin[tls](4),prime(bin[tls](3)),bin[l_past](phot),prime(bin[l_past](phot+1)),-gam_r_fb*sqrt(1.*phot));        
    // the adjoint left atom interaction hamiltonian
    H_fb_r.set(bin[tls](1),prime(bin[tls](2)),bin[l_past](phot+1),prime(bin[l_past](phot)), c_gam_r_fb*sqrt(1.*phot));
    H_fb_r.set(bin[tls](1),prime(bin[tls](3)),bin[l_past](phot+1),prime(bin[l_past](phot)), c_gam_r_fb*sqrt(1.*phot));        
    H_fb_r.set(bin[tls](2),prime(bin[tls](4)),bin[l_past](phot+1),prime(bin[l_past](phot)), c_gam_r_fb*sqrt(1.*phot));        
    H_fb_r.set(bin[tls](3),prime(bin[tls](4)),bin[l_past](phot+1),prime(bin[l_past](phot)),-c_gam_r_fb*sqrt(1.*phot));        
    }
    auto H_fb  =   H_fb_l * delta(bin[l_past],prime(bin[l_past])) 
                 + H_fb_r * delta(bin[r_past],prime(bin[r_past]))
               ;

    auto H_int =  delta(bin[l_past],prime(bin[l_past])) * H_dis * delta(bin[r_past],prime(bin[r_past])) 
                + delta(bin[l_now ],prime(bin[l_now ])) * H_fb  * delta(bin[r_now ],prime(bin[r_now ]))    
               ;
    //PrintData(H_int);           

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
 
    auto temp= ITensor(bin[tls],prime(bin[tls]));
    temp.set(bin[tls](1),prime(bin[tls](1)),1.);
    temp.set(bin[tls](2),prime(bin[tls](2)),1.);
    temp.set(bin[tls](3),prime(bin[tls](3)),1.);
    temp.set(bin[tls](4),prime(bin[tls](4)),1.);

    temp = delta(bin[l_now],prime(bin[l_now]))*temp*delta(bin[l_past],prime(bin[l_past]));
    temp = delta(bin[r_now],prime(bin[r_now]))*temp*delta(bin[r_past],prime(bin[r_past]));
    // now the taylor expansion of the U_evo
    auto U_evo = temp + H_int_1 + H_int_2 + H_int_3 + H_int_4 + H_int_5 + H_int_6 + H_int_7  + H_int_8 + H_int_9 + H_int_10 ;    
    //PrintData(U_evo);
    
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // --------------------- END U_EVO SETUP --------------------------------------------
    // ----------------------------------------------------------------------------------
    double pop1=0.;
    double pop2=0.;
    double mps_norm=1.;
    
    ITensor U,S,V,W,SWAP;
    
// das bild im MPS sieht aus [pastbin][feedbackbin][systembin][currentbin]
Index iFBl,iFBr,iCBl,iCBr; // index fuer feedback, system, und current bin 
Index iLinkp,iLinkc; // index fuer feedback bin und past bin

int l; 
//for(l=1;l<=ttotal+1;l++) {printf("%i:\n",l); PrintData(psi(l));}
// psi(1) == bin[1], ... ,psi(Nfb) == bin[Nfb], psi(Nfb+1) == bin[0], psi(Nfb+2) == bin[Nfb+1], .... psi(ttotal+1) = bin[ttotal] 
// sys_at=Nfb+1 at the beginning
// algorithm moves bin[0] forward, and the feedback distance is Nfb=2*fb
    
    time_t curtime; struct tm *loctime; curtime = time(NULL); loctime = localtime (&curtime); 
    int hour = loctime -> tm_hour;int minute = loctime -> tm_min;int second = loctime -> tm_sec;
    int year = 1900+loctime -> tm_year;int month = loctime -> tm_mon;int day  = loctime -> tm_mday;
    printf("Date: %d.%d.%.d -- Time: %d:%d:%d  \n",day,month+1,year,hour,minute,second);  

    FILE * file;
    char FILE_NAME[2048+1024];
    snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_TLS_w_FEEDBACK_at_Step=%i_of%i_init1_%.2f_init2_%.2f_init3_%.2f_init4_%.2f_Gam_l%.3f_Gam_r%.3f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,fb,t_end,init_GG,init_SYM,init_ASYM,init_EE,Gamma1,Gamma2,dt,cutoff*1000000.);
    file = fopen(FILE_NAME,"w");
    
    pop1= ( state2_pop(psi(sys_at))+state3_pop(psi(sys_at))+2.*state4_pop(psi(sys_at))-2.*re_coherence23(psi(sys_at)) )*0.5;
    pop2= ( state2_pop(psi(sys_at))+state3_pop(psi(sys_at))+2.*state4_pop(psi(sys_at))+2.*re_coherence23(psi(sys_at)) )*0.5;
    mps_norm = TLS_norm(psi(sys_at));

for(int m=0;m<t_end;m++)
{   
    printf("Step %i: pop1=%.10f -- pop2=%.10f -- norm=%.10f -- OrthoCenter=%i \n",m+1,pop1,pop2,mps_norm,0);
    fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \n",(m+1.)*dt,pop1,pop2,mps_norm);

    // --------------------------------------------------------
    // ----------- SWAP TO U_EVO START ------------------------
    for(int k=0;k<Nfb-2;k++)
    {    
    l_past=sys_at-Nfb  +k;
    r_past=sys_at-Nfb+1+k;
    
    SWAP = psi(r_past)*psi(r_past+1);
    U=ITensor(bin[r_past+1],commonIndex(psi(r_past-1),psi(r_past)));
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});
    psi.setA(r_past,U); 
    psi.setA(r_past+1,S*V);
    
    SWAP = psi(l_past)*psi(r_past);
    if ( l_past>1 ) U=ITensor(bin[r_past+1],commonIndex(psi(l_past),psi(l_past-1)));
    else U=ITensor(bin[r_past+1]);
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});
    psi.setA(l_past,U); 
    psi.setA(l_past+1,S*V);
    }
    // ----------- SWAP TO U_EVO END --------------------------
    // --------------------------------------------------------

    // --------------------------------------------------------
    // ----------- APPLY MPO START ----------------------------
    iFBl = findIndex(psi.A(sys_at-2),"left"); 
    iFBr = findIndex(psi.A(sys_at-1),"right");      
    iCBl = findIndex(psi.A(sys_at+1),"left");
    iCBr = findIndex(psi.A(sys_at+2),"right");    
    iLinkp = commonIndex(psi.A(sys_at-2),psi.A(sys_at-3));
    iLinkc = commonIndex(psi.A(sys_at+2),psi.A(sys_at+3));
    
    temp = noPrime(psi(sys_at-2)*psi(sys_at-1)*U_evo*psi(sys_at)*psi(sys_at+1)*psi(sys_at+2)); //PrintData(temp);

    U=ITensor(iLinkp,iFBl,iFBr,iCBl,iCBr); 
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    // tls in the next mps slot 
    psi.setA(sys_at+2,V); //psi.position(sys_at+2);
 
    pop1= ( state2_pop(V*S)+state3_pop(V*S)+2.*state4_pop(V*S)-2.*re_coherence23(V*S) )*0.5;
    pop2= ( state2_pop(V*S)+state3_pop(V*S)+2.*state4_pop(V*S)+2.*re_coherence23(V*S) )*0.5;
    mps_norm = TLS_norm(V*S);
    
    // now factorize step by step the reservoir tensor
    temp=U*S; 
    U=ITensor(iLinkp,iFBl,iFBr,iCBl); 
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    psi.setA(sys_at+1,V); // iCBr 
    
    temp=U*S; 
    U=ITensor(iLinkp,iFBl,iFBr); 
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    psi.setA(sys_at,V); // iCBl 

    temp=U*S; 
    U=ITensor(iLinkp,iFBl); 
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    psi.setA(sys_at-1,V); // iFBr 
    psi.setA(sys_at-2,U*S); // iFBl
    
    //PrintData(psi(sys_at-2));    PrintData(psi(sys_at-1));     PrintData(psi(sys_at));
    //PrintData(psi(sys_at+1));    PrintData(psi(sys_at+2));
    // ----------- APPLY MPO END  -----------------------------
    // --------------------------------------------------------
      
    // --------------------------------------------------------
    // ----------- SWAP BACK START ----------------------------
    for(int k=0; k<Nfb-2;k++)
    {    
    l_past=sys_at-2-k;
    r_past=sys_at-1-k;

    SWAP = psi(l_past)*psi(l_past-1);
    U=ITensor(bin[sys_at-1-k],commonIndex(psi(l_past),psi(r_past)));
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});
    if (k==Nfb-3) { psi.setA(l_past,U*S); psi.setA(l_past-1,V); } // if -- else necessary to store
    else          { psi.setA(l_past,U); psi.setA(l_past-1,S*V); } // the normalization in the next pastbin

    SWAP = psi(r_past)*psi(l_past);
    U=ITensor(bin[sys_at-1-k],commonIndex(psi(r_past),psi(r_past+1)));  
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});              
    if (k==Nfb-3) { psi.setA(r_past,U*S);  psi.setA(l_past,V); } // if -- else necessary to store
    else          { psi.setA(r_past,U);  psi.setA(l_past,S*V); } // the normalization in the next pastbin
    }
    
    // ----------- SWAP BACK END ------------------------------
    // --------------------------------------------------------
    // Update MPO, only if a further loop is possible
    if (m<t_end-1)
    {    
    l_past = 1 + m*2; r_past = 2 + m*2; l_now = Nfb+1+m*2; r_now = Nfb+2+m*2;
    U_evo = delta(bin[l_past],bin[l_past+2])*U_evo*delta(prime(bin[l_past]),prime(bin[l_past+2])) ;
    U_evo = delta(bin[r_past],bin[r_past+2])*U_evo*delta(prime(bin[r_past]),prime(bin[r_past+2])) ;
    U_evo = delta(bin[l_now],bin[l_now+2])*U_evo*delta(prime(bin[l_now]),prime(bin[l_now+2])) ;
    U_evo = delta(bin[r_now],bin[r_now+2])*U_evo*delta(prime(bin[r_now]),prime(bin[r_now+2])) ;
    //PrintData(U_evo);
    }
    // increase the position of the system by 2!!
    sys_at++;sys_at++;
}   
//printf("Here!\n"); //    for(l=1;l<=ttotal+1;l++) {printf("%i:\n",l); PrintData(psi(l));}    

int ehour = loctime -> tm_hour; int eminute = loctime -> tm_min; int esecond = loctime -> tm_sec;
int eyear = 1900+loctime -> tm_year; int emonth = loctime -> tm_mon; int eday  = loctime -> tm_mday;
double time_before,time_now; time_before =esecond + 60*eminute+60*60*ehour;

// give information how long it took
curtime = time(NULL); loctime = localtime (&curtime); 
ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
eyear = 1900+loctime -> tm_year; emonth = loctime -> tm_mon; eday  = loctime -> tm_mday;
printf("End of calculation: %d.%d.%.d -- Time: %d:%d:%d  \n",eday,emonth+1,eyear,ehour,eminute,esecond);  
fprintf(file,"# end of calculation %02d_%02d_%02d_%02d_%02d_%02d \n",eyear,emonth+1,eday,ehour,eminute,esecond);
printf("Time elapsed: %d.%d.%.d -- Time: %dh:%dmin:%dsec \n",eday-day,emonth-month,eyear-year,ehour-hour,eminute-minute,esecond-second);  
fprintf(file,"# time difference %02d_%02d_%02d_%02d_%02d_%02d \n",eyear-year,emonth-month,eday-day,ehour-hour,eminute-minute,esecond-second);


fclose(file);

return 0;
}
