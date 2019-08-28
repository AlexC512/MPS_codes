#include "itensor/all.h"

/*

Below is a very inefficient code, because I have not 
the standard form of an mps where every bin has one physical index and two link indices,
this allows fast svd and fast swap protocols. example one and same calculation with 

time_end = 35;
Gamma1 =0.20;
Gamma2 =0.20;
Nbin =2;
svdcutoff = 1E-8;
maxnumberofSV = 5000;
feedback_time = 20;
init_EE   = 1.0;

takes with a standard mps 14 seconds
with the code below 1 minute 11 seconds

*/

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
    int Nfb = fb;
    // size of array is twice because of left and right moving photons and one feedback loop 
    int ttotal = t_end+Nfb; printf("Total Number of time bins: %i\n",ttotal);
    
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
    auto lbin = std::vector<Index>((int)ttotal+1); 
    auto rbin = std::vector<Index>((int)ttotal+1); 
    for(int j = 1; j <= ttotal; ++j){ lbin.at(j) = Index(3,"left");     }
    for(int j = 1; j <= ttotal; ++j){ rbin.at(j) = Index(3,"right");     }
    int tls=0;
    lbin.at(tls) = Index(4,"sys"); // bin(0) system consists of two atoms 
        
    // increase number by one because we do not use zeroth link!!
    auto binlink = std::vector<Index>((int)ttotal+1);
    for(int j = 1; j <= ttotal; ++j) { binlink.at(j) = Index(1,"Link");  }
    // we need space in the mps ttotal time bins + 1 system bin ... actually there is psi(0)
    // but it does not count?
    MPS psi=MPS(ttotal+1);
    // set up the mps
    // I need fb bins before system and current bin
    // psi(1) ... psi(fb) ... then psi(fb+1)=sys, psi(fb+2) current bin
    // bin[1] ... bin[fb] ... then bin[0], bin[fb+1], bin[fb+2] ... bin[fb+t_end] 
    // now=fb+1, past=1=now-fb, to remember psi(now) is the system, bin[now] is in psi(now+1)
     ITensor MPSTensor = ITensor();
     // first MPS slot is filled with a bin which has no link to the left
     MPSTensor = ITensor(lbin[1],rbin[1],binlink[1]);
     MPSTensor.set(lbin[1](1),rbin[1](1),binlink[1](1),1.); // Level besetzt im TLS
     psi.setA((1),MPSTensor);

    for(int j = 2; j<=fb; j++)
    {
     MPSTensor = ITensor(lbin[j],rbin[j],binlink[j-1],binlink[j]);
     MPSTensor.set(lbin[j](1),rbin[j](1),binlink[j-1](1),binlink[j](1),1.); // Level besetzt im TLS
     psi.setA((j),MPSTensor);
    }   

     MPSTensor = ITensor(lbin[tls],binlink[fb],binlink[fb+1]);
     MPSTensor.set(lbin[tls](1),binlink[fb](1),binlink[fb+1](1),init_GG);   
     MPSTensor.set(lbin[tls](2),binlink[fb](1),binlink[fb+1](1),init_SYM);   
     MPSTensor.set(lbin[tls](3),binlink[fb](1),binlink[fb+1](1),init_ASYM);   
     MPSTensor.set(lbin[tls](4),binlink[fb](1),binlink[fb+1](1),init_EE);   
     psi.setA((fb+1),MPSTensor);

     for(int j = fb+1; j<ttotal; j++)
    {
     MPSTensor = ITensor(lbin[j],rbin[j],binlink[j],binlink[j+1]);
     MPSTensor.set(lbin[j](1),rbin[j](1),binlink[j](1),binlink[j+1](1),1.); // Level besetzt im TLS
     psi.setA((j+1),MPSTensor);
    }   

    MPSTensor = ITensor(lbin[ttotal],rbin[ttotal],binlink[ttotal]);
    MPSTensor.set(lbin[ttotal](1),rbin[ttotal](1),binlink[ttotal](1),1.); // Level besetzt im TLS
    psi.setA((ttotal+1),MPSTensor);
    // now we have created an mps with (2*t_end time+Nfb) bins + 1 system bin located so
    // that feedback bins are initially empty to the left
    
    psi.position(1);
    printf("Orthocenter at psi(%i) \n",orthoCenter(psi));

    int l; //    for(l=1;l<=ttotal+1;l++) {printf("%i:\n",l); PrintData(psi(l));}
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ---------------------END MPS SETUP -----------------------------------------------
    // ----------------------------------------------------------------------------------
  
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ---------------------U_EVO SETUP -------------------------------------------------
    // ----------------------------------------------------------------------------------
    int now = fb+1;
        
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
   
    auto H_dis_l = ITensor(lbin[tls],prime(lbin[tls]),lbin[now],prime(lbin[now]));
    auto H_dis_r = ITensor(lbin[tls],prime(lbin[tls]),rbin[now],prime(rbin[now]));
    for( int phot=1;phot<3;phot++)
    {    
    H_dis_l.set(lbin[tls](2),prime(lbin[tls](1)),lbin[now](phot),prime(lbin[now](phot+1)), gam_l*sqrt(1.*phot));
    H_dis_l.set(lbin[tls](3),prime(lbin[tls](1)),lbin[now](phot),prime(lbin[now](phot+1)),-gam_l*sqrt(1.*phot));        
    H_dis_l.set(lbin[tls](4),prime(lbin[tls](2)),lbin[now](phot),prime(lbin[now](phot+1)), gam_l*sqrt(1.*phot));        
    H_dis_l.set(lbin[tls](4),prime(lbin[tls](3)),lbin[now](phot),prime(lbin[now](phot+1)), gam_l*sqrt(1.*phot));        
    // the adjoint left atom interaction hamiltonian
    H_dis_l.set(lbin[tls](1),prime(lbin[tls](2)),lbin[now](phot+1),prime(lbin[now](phot)), c_gam_l*sqrt(1.*phot));
    H_dis_l.set(lbin[tls](1),prime(lbin[tls](3)),lbin[now](phot+1),prime(lbin[now](phot)),-c_gam_l*sqrt(1.*phot));        
    H_dis_l.set(lbin[tls](2),prime(lbin[tls](4)),lbin[now](phot+1),prime(lbin[now](phot)), c_gam_l*sqrt(1.*phot));        
    H_dis_l.set(lbin[tls](3),prime(lbin[tls](4)),lbin[now](phot+1),prime(lbin[now](phot)), c_gam_l*sqrt(1.*phot));        
    // right atom interation hamiltonian ... photons in time bins are created
    H_dis_r.set(lbin[tls](2),prime(lbin[tls](1)),rbin[now](phot),prime(rbin[now](phot+1)), gam_r*sqrt(1.*phot));
    H_dis_r.set(lbin[tls](3),prime(lbin[tls](1)),rbin[now](phot),prime(rbin[now](phot+1)), gam_r*sqrt(1.*phot));        
    H_dis_r.set(lbin[tls](4),prime(lbin[tls](2)),rbin[now](phot),prime(rbin[now](phot+1)), gam_r*sqrt(1.*phot));        
    H_dis_r.set(lbin[tls](4),prime(lbin[tls](3)),rbin[now](phot),prime(rbin[now](phot+1)),-gam_r*sqrt(1.*phot));        
    // the adjoint right atom interaction hamiltonian ... photons in time bins are annihilated
    H_dis_r.set(lbin[tls](1),prime(lbin[tls](2)),rbin[now](phot+1),prime(rbin[now](phot)), c_gam_r*sqrt(1.*phot));
    H_dis_r.set(lbin[tls](1),prime(lbin[tls](3)),rbin[now](phot+1),prime(rbin[now](phot)), c_gam_r*sqrt(1.*phot));        
    H_dis_r.set(lbin[tls](2),prime(lbin[tls](4)),rbin[now](phot+1),prime(rbin[now](phot)), c_gam_r*sqrt(1.*phot));        
    H_dis_r.set(lbin[tls](3),prime(lbin[tls](4)),rbin[now](phot+1),prime(rbin[now](phot)),-c_gam_r*sqrt(1.*phot));        
    }
    auto H_dis =   H_dis_l * delta(rbin[now],prime(rbin[now])) 
                 + H_dis_r * delta(lbin[now],prime(lbin[now])) 
               ;  

   // time non-local field with excitation exchange ... ATTENTION ... single excitation limit 
    // interaction with the pastbin
    auto H_fb_l  = ITensor(lbin[tls],prime(lbin[tls]),rbin[now-fb],prime(rbin[now-fb]));
    auto H_fb_r  = ITensor(lbin[tls],prime(lbin[tls]),lbin[now-fb],prime(lbin[now-fb]));
    for( int phot=1;phot<3;phot++)
    {    
    H_fb_l.set(lbin[tls](2),prime(lbin[tls](1)),rbin[now-fb](phot),prime(rbin[now-fb](phot+1)), gam_l_fb*sqrt(1.*phot));
    H_fb_l.set(lbin[tls](3),prime(lbin[tls](1)),rbin[now-fb](phot),prime(rbin[now-fb](phot+1)),-gam_l_fb*sqrt(1.*phot));        
    H_fb_l.set(lbin[tls](4),prime(lbin[tls](2)),rbin[now-fb](phot),prime(rbin[now-fb](phot+1)), gam_l_fb*sqrt(1.*phot));        
    H_fb_l.set(lbin[tls](4),prime(lbin[tls](3)),rbin[now-fb](phot),prime(rbin[now-fb](phot+1)), gam_l_fb);        
    // the adjoint left atom interaction hamiltonian
    H_fb_l.set(lbin[tls](1),prime(lbin[tls](2)),rbin[now-fb](phot+1),prime(rbin[now-fb](phot)), c_gam_l_fb*sqrt(1.*phot));
    H_fb_l.set(lbin[tls](1),prime(lbin[tls](3)),rbin[now-fb](phot+1),prime(rbin[now-fb](phot)),-c_gam_l_fb*sqrt(1.*phot));        
    H_fb_l.set(lbin[tls](2),prime(lbin[tls](4)),rbin[now-fb](phot+1),prime(rbin[now-fb](phot)), c_gam_l_fb*sqrt(1.*phot));        
    H_fb_l.set(lbin[tls](3),prime(lbin[tls](4)),rbin[now-fb](phot+1),prime(rbin[now-fb](phot)), c_gam_l_fb*sqrt(1.*phot));                
    // time non-local field ... ATTENTION ... single excitation limit 
    H_fb_r.set(lbin[tls](2),prime(lbin[tls](1)),lbin[now-fb](phot),prime(lbin[now-fb](phot+1)), gam_r_fb*sqrt(1.*phot));
    H_fb_r.set(lbin[tls](3),prime(lbin[tls](1)),lbin[now-fb](phot),prime(lbin[now-fb](phot+1)), gam_r_fb*sqrt(1.*phot));        
    H_fb_r.set(lbin[tls](4),prime(lbin[tls](2)),lbin[now-fb](phot),prime(lbin[now-fb](phot+1)), gam_r_fb*sqrt(1.*phot));        
    H_fb_r.set(lbin[tls](4),prime(lbin[tls](3)),lbin[now-fb](phot),prime(lbin[now-fb](phot+1)),-gam_r_fb*sqrt(1.*phot));        
    // the adjoint left atom interaction hamiltonian
    H_fb_r.set(lbin[tls](1),prime(lbin[tls](2)),lbin[now-fb](phot+1),prime(lbin[now-fb](phot)), c_gam_r_fb*sqrt(1.*phot));
    H_fb_r.set(lbin[tls](1),prime(lbin[tls](3)),lbin[now-fb](phot+1),prime(lbin[now-fb](phot)), c_gam_r_fb*sqrt(1.*phot));        
    H_fb_r.set(lbin[tls](2),prime(lbin[tls](4)),lbin[now-fb](phot+1),prime(lbin[now-fb](phot)), c_gam_r_fb*sqrt(1.*phot));        
    H_fb_r.set(lbin[tls](3),prime(lbin[tls](4)),lbin[now-fb](phot+1),prime(lbin[now-fb](phot)),-c_gam_r_fb*sqrt(1.*phot));        
    }
    auto H_fb  =   H_fb_l * delta(lbin[now-fb],prime(lbin[now-fb])) 
                 + H_fb_r * delta(rbin[now-fb],prime(rbin[now-fb]))
               ;

    auto H_int =  delta(lbin[now-fb],prime(lbin[now-fb])) * H_dis * delta(rbin[now-fb],prime(rbin[now-fb])) 
                + delta(lbin[now ],prime(lbin[now ])) * H_fb  * delta(rbin[now ],prime(rbin[now ]))    
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

    auto temp= ITensor(lbin[tls],prime(lbin[tls]));
    temp.set(lbin[tls](1),prime(lbin[tls](1)),1.);
    temp.set(lbin[tls](2),prime(lbin[tls](2)),1.);
    temp.set(lbin[tls](3),prime(lbin[tls](3)),1.);
    temp.set(lbin[tls](4),prime(lbin[tls](4)),1.);

    temp = delta(lbin[now],prime(lbin[now]))*temp*delta(lbin[now-fb],prime(lbin[now-fb]));
    temp = delta(rbin[now],prime(rbin[now]))*temp*delta(rbin[now-fb],prime(rbin[now-fb]));
    // now the taylor expansion of the U_evo
    auto U_evo = temp 
                 + H_int_1 + H_int_2 + H_int_3 + H_int_4 + H_int_5 + H_int_6 + H_int_7  + H_int_8 + H_int_9 + H_int_10 ;    
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

    time_t curtime; struct tm *loctime; curtime = time(NULL); loctime = localtime (&curtime); 
    int hour = loctime -> tm_hour;int minute = loctime -> tm_min;int second = loctime -> tm_sec;
    int year = 1900+loctime -> tm_year;int month = loctime -> tm_mon;int day  = loctime -> tm_mday;
    printf("Date: %d.%d.%.d -- Time: %d:%d:%d  \n",day,month+1,year,hour,minute,second);  

    FILE * file;
    char FILE_NAME[2048+1024];
    snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_TLS_w_FEEDBACK_at_Step=%i_of%i_init1_%.2f_init2_%.2f_init3_%.2f_init4_%.2f_Gam_l%.3f_Gam_r%.3f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,fb,t_end,init_GG,init_SYM,init_ASYM,init_EE,Gamma1,Gamma2,dt,cutoff*1000000.);
    file = fopen(FILE_NAME,"w");
    
    pop1= ( state2_pop(psi(now))+state3_pop(psi(now))+2.*state4_pop(psi(now))-2.*re_coherence23(psi(now)) )*0.5;
    pop2= ( state2_pop(psi(now))+state3_pop(psi(now))+2.*state4_pop(psi(now))+2.*re_coherence23(psi(now)) )*0.5;
    mps_norm = TLS_norm(psi(now));

    // psi(1) ... psi(fb) psi(fb+1) psi(fb+2) ... psi(fb+t_end)
    // you swap psi(1)=psi(now-fb) to psi(now-1), then contract psi(now-1)*psi(now)
    // the tensor is now ready for the mpo, apply mpo, 
    //create psi(now-1) with pastbin, psi(now) with currentbin, psi(now+1) with systembin 
    // dynamics psi(1) psi(2) ... psi(fb+1) psi( 
    
for(int m=0;m<t_end;m++)
{   
    printf("Step %i: pop1=%.10f -- pop2=%.10f -- norm=%.10f \n",m+1,pop1,pop2,mps_norm);
    fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \n",(m+1.)*dt,pop1,pop2,mps_norm);

    // --------------------------------------------------------
    // ----------- SWAP TO U_EVO START ------------------------
    for(int k=now-fb;k<now-1;k++)
    {        
    SWAP = psi(k)*psi(k+1);
    if (k==1)  U=ITensor(lbin[k+1],rbin[k+1]);
    else       U=ITensor(lbin[k+1],rbin[k+1],commonIndex(psi(k-1),psi(k)));
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});
    psi.setA(k,U); 
    psi.setA(k+1,S*V);
    }
    // ----------- SWAP TO U_EVO END --------------------------
    // --------------------------------------------------------

    
    // --------------------------------------------------------
    // ----------- APPLY MPO START ----------------------------
    iFBl = findIndex(psi.A(now-1),"left"); 
    iFBr = findIndex(psi.A(now-1),"right");      
    iCBl = findIndex(psi.A(now+1),"left");
    iCBr = findIndex(psi.A(now+1),"right");    
    iLinkp = commonIndex(psi.A(now-1),psi.A(now-2));
    iLinkc = commonIndex(psi.A(now+1),psi.A(now+2));
    
    temp = noPrime(psi(now-1)*U_evo*psi(now)*psi(now+1)); 

    U=ITensor(iLinkp,iFBl,iFBr,iCBl,iCBr); 
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    psi.setA(now+1,V);
 
    pop1= ( state2_pop(V*S)+state3_pop(V*S)+2.*state4_pop(V*S)-2.*re_coherence23(V*S) )*0.5;
    pop2= ( state2_pop(V*S)+state3_pop(V*S)+2.*state4_pop(V*S)+2.*re_coherence23(V*S) )*0.5;
    mps_norm = TLS_norm(V*S);
    
    // now factorize step by step the reservoir tensor
    temp=U*S; 
    U=ITensor(iLinkp,iFBl,iFBr); 
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    psi.setA(now,V); 
    psi.setA(now-1,U*S);    
    // ----------- APPLY MPO END  -----------------------------
    // --------------------------------------------------------
      
    // --------------------------------------------------------
    // ----------- SWAP BACK START ----------------------------    
    for(int k=now-1;k>now-fb;k--)
    {        
    SWAP = psi(k)*psi(k-1); 
    U=ITensor(lbin[now-fb],rbin[now-fb],commonIndex(psi(k-1),psi(k-2)));
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});
    if (k==now-fb-1) { psi.setA(k-1,U*S); psi.setA(k,V); }
    else             { psi.setA(k-1,U); psi.setA(k,V*S); }
    }
    // ----------- SWAP BACK END ------------------------------
    // --------------------------------------------------------
    // Update MPO, only if a further loop is possible
    if (m<t_end-1)
    {    
    U_evo = delta(lbin[now-fb],lbin[now-fb+1])*U_evo*delta(prime(lbin[now-fb]),prime(lbin[now-fb+1])) ;
    U_evo = delta(rbin[now-fb],rbin[now-fb+1])*U_evo*delta(prime(rbin[now-fb]),prime(rbin[now-fb+1])) ;
    U_evo = delta(lbin[now],lbin[now+1])*U_evo*delta(prime(lbin[now]),prime(lbin[now+1])) ;
    U_evo = delta(rbin[now],rbin[now+1])*U_evo*delta(prime(rbin[now]),prime(rbin[now+1])) ;
    }
    now++;
    
}   

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
