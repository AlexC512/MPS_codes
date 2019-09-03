#include "itensor/all.h"

using namespace itensor;


void MPO_SETUP(ITensor& U_evo, const std::vector<Index>& bin,
               int tls, int l_past, int l_middle, int l_now, int Nbin,
               Real Gamma_l, Real Gamma_m, Real Gamma_r,
               Real phi_l, Real phi_m, Real phi_r, Real dt
              )
{
    
    //int l_past= sys_at-Nfb;
    int r_past= l_past+1;
    //int l_now = Nfb+1; 
    int r_now = l_now+1;
    //int l_middle = sys_at-Nfb/2   ;
    int r_middle = l_middle +1;
    
    // left atom interaction hamiltonian ... time local field (now)
    // define: |abc> == Atom left in a, Atom middle b, Atom right c, Atom 3, Atom 2, Atom 1, 
    Complex   gam_l    =  -Cplx_i*Gamma_l*sqrt(dt);
    Complex c_gam_l    =  -Cplx_i*Gamma_l*sqrt(dt);
    Complex   gam_l_fb =  -Cplx_i*Gamma_l*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_l);
    Complex c_gam_l_fb =  -Cplx_i*Gamma_l*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_l);

    Complex   gam_m_fb =  -Cplx_i*Gamma_m*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_m); // ATTENTION change phase definition
    Complex c_gam_m_fb =  -Cplx_i*Gamma_m*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_m); // if left and right feedback is not the same!!!

    Complex   gam_r    =  -Cplx_i*Gamma_r*sqrt(dt);
    Complex c_gam_r    =  -Cplx_i*Gamma_r*sqrt(dt);
    Complex   gam_r_fb =  -Cplx_i*Gamma_r*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_r);
    Complex c_gam_r_fb =  -Cplx_i*Gamma_r*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_r);
   
    auto H_dis_l  = ITensor(bin[tls],prime(bin[tls]),bin[l_now],prime(bin[l_now]));
    auto H_dis_ml = ITensor(bin[tls],prime(bin[tls]),bin[l_middle],prime(bin[l_middle]));
    auto H_dis_mr = ITensor(bin[tls],prime(bin[tls]),bin[r_middle],prime(bin[r_middle]));
    auto H_dis_r  = ITensor(bin[tls],prime(bin[tls]),bin[r_now],prime(bin[r_now]));
    double sqrt_phot=1.;
    for( int phot=1;phot<Nbin;phot++)
    {    
    sqrt_phot=sqrt(1.*phot);    
    // -------------------------------- left atom -----------------------------------------------------------------
    // for the left atom deexcitation means -4 we go from unprimed to primed    
    H_dis_l.set(bin[tls](8),prime(bin[tls](4)),bin[l_now](phot),prime(bin[l_now](phot+1)),gam_l*sqrt_phot);
    H_dis_l.set(bin[tls](7),prime(bin[tls](3)),bin[l_now](phot),prime(bin[l_now](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](6),prime(bin[tls](2)),bin[l_now](phot),prime(bin[l_now](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](5),prime(bin[tls](1)),bin[l_now](phot),prime(bin[l_now](phot+1)),gam_l*sqrt_phot);        
    // for the left atom excitation means +4
    H_dis_l.set(bin[tls](4),prime(bin[tls](8)),bin[l_now](phot+1),prime(bin[l_now](phot)), c_gam_l*sqrt_phot);
    H_dis_l.set(bin[tls](3),prime(bin[tls](7)),bin[l_now](phot+1),prime(bin[l_now](phot)), c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](2),prime(bin[tls](6)),bin[l_now](phot+1),prime(bin[l_now](phot)), c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](1),prime(bin[tls](5)),bin[l_now](phot+1),prime(bin[l_now](phot)), c_gam_l*sqrt_phot);        
    // -------------------------------- middle atom -----------------------------------------------------------------
    // middle atom interation hamiltonian ... photons in time bins are created ... middle atom has no separate feedback int
    H_dis_ml.set(bin[tls](8),prime(bin[tls](6)),bin[l_middle](phot),prime(bin[l_middle](phot+1)), gam_m_fb*sqrt_phot);
    H_dis_ml.set(bin[tls](7),prime(bin[tls](5)),bin[l_middle](phot),prime(bin[l_middle](phot+1)), gam_m_fb*sqrt_phot);        
    H_dis_ml.set(bin[tls](4),prime(bin[tls](2)),bin[l_middle](phot),prime(bin[l_middle](phot+1)), gam_m_fb*sqrt_phot);        
    H_dis_ml.set(bin[tls](3),prime(bin[tls](1)),bin[l_middle](phot),prime(bin[l_middle](phot+1)), gam_m_fb*sqrt_phot);        
    // interaction with left moving field
    H_dis_ml.set(bin[tls](6),prime(bin[tls](8)),bin[l_middle](phot+1),prime(bin[l_middle](phot)), c_gam_m_fb*sqrt_phot);
    H_dis_ml.set(bin[tls](5),prime(bin[tls](7)),bin[l_middle](phot+1),prime(bin[l_middle](phot)), c_gam_m_fb*sqrt_phot);        
    H_dis_ml.set(bin[tls](2),prime(bin[tls](4)),bin[l_middle](phot+1),prime(bin[l_middle](phot)), c_gam_m_fb*sqrt_phot);        
    H_dis_ml.set(bin[tls](1),prime(bin[tls](3)),bin[l_middle](phot+1),prime(bin[l_middle](phot)), c_gam_m_fb*sqrt_phot); 
    // middle atom interation hamiltonian ... photons in time bins are created ... setting the location x zero
    H_dis_mr.set(bin[tls](8),prime(bin[tls](6)),bin[r_middle](phot),prime(bin[r_middle](phot+1)), gam_m_fb*sqrt_phot);
    H_dis_mr.set(bin[tls](7),prime(bin[tls](5)),bin[r_middle](phot),prime(bin[r_middle](phot+1)), gam_m_fb*sqrt_phot);        
    H_dis_mr.set(bin[tls](4),prime(bin[tls](2)),bin[r_middle](phot),prime(bin[r_middle](phot+1)), gam_m_fb*sqrt_phot);        
    H_dis_mr.set(bin[tls](3),prime(bin[tls](1)),bin[r_middle](phot),prime(bin[r_middle](phot+1)), gam_m_fb*sqrt_phot);        
    // interaction with right moving field
    H_dis_mr.set(bin[tls](6),prime(bin[tls](8)),bin[r_middle](phot+1),prime(bin[r_middle](phot)), c_gam_m_fb*sqrt_phot);
    H_dis_mr.set(bin[tls](5),prime(bin[tls](7)),bin[r_middle](phot+1),prime(bin[r_middle](phot)), c_gam_m_fb*sqrt_phot);        
    H_dis_mr.set(bin[tls](2),prime(bin[tls](4)),bin[r_middle](phot+1),prime(bin[r_middle](phot)), c_gam_m_fb*sqrt_phot);        
    H_dis_mr.set(bin[tls](1),prime(bin[tls](3)),bin[r_middle](phot+1),prime(bin[r_middle](phot)), c_gam_m_fb*sqrt_phot); 
    // -------------------------------- right atom -----------------------------------------------------------------
    // for the left atom deexcitation means -4 we go from unprimed to primed    
    H_dis_r.set(bin[tls](8),prime(bin[tls](7)),bin[r_now](phot),prime(bin[r_now](phot+1)), gam_r*sqrt_phot);
    H_dis_r.set(bin[tls](6),prime(bin[tls](5)),bin[r_now](phot),prime(bin[r_now](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](4),prime(bin[tls](3)),bin[r_now](phot),prime(bin[r_now](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](2),prime(bin[tls](1)),bin[r_now](phot),prime(bin[r_now](phot+1)), gam_r*sqrt_phot);        
    // for the left atom excitation means +4
    H_dis_r.set(bin[tls](7),prime(bin[tls](8)),bin[r_now](phot+1),prime(bin[r_now](phot)), c_gam_r*sqrt_phot);
    H_dis_r.set(bin[tls](5),prime(bin[tls](6)),bin[r_now](phot+1),prime(bin[r_now](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](3),prime(bin[tls](4)),bin[r_now](phot+1),prime(bin[r_now](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](1),prime(bin[tls](2)),bin[r_now](phot+1),prime(bin[r_now](phot)), c_gam_r*sqrt_phot);        
    }
    
    // now bring dissipation hamiltonians into one product Hilbert space
    H_dis_l = delta(bin[r_now],prime(bin[r_now]))       * H_dis_l;
    H_dis_l = delta(bin[r_middle],prime(bin[r_middle])) * H_dis_l;
    H_dis_l = delta(bin[l_middle],prime(bin[l_middle])) * H_dis_l;
    
    H_dis_ml = delta(bin[r_now],prime(bin[r_now]))       * H_dis_ml;
    H_dis_ml = delta(bin[r_middle],prime(bin[r_middle])) * H_dis_ml;
    H_dis_ml = delta(bin[l_now],prime(bin[l_now]))       * H_dis_ml;
    
    H_dis_mr = delta(bin[r_now],prime(bin[r_now]))       * H_dis_mr;
    H_dis_mr = delta(bin[l_middle],prime(bin[l_middle])) * H_dis_mr;
    H_dis_mr = delta(bin[l_now],prime(bin[l_now]))       * H_dis_mr;

    H_dis_r = delta(bin[l_now],prime(bin[l_now]))       * H_dis_r;
    H_dis_r = delta(bin[r_middle],prime(bin[r_middle])) * H_dis_r;
    H_dis_r = delta(bin[l_middle],prime(bin[l_middle])) * H_dis_r;
    
    auto H_dis =   H_dis_l + H_dis_ml + H_dis_mr + H_dis_r    ;  
   
    // time non-local field with excitation exchange for left and right atom as the middle one has no extra feedback 
    auto H_fb_l  = ITensor(bin[tls],prime(bin[tls]),bin[r_past],prime(bin[r_past]));
    auto H_fb_r  = ITensor(bin[tls],prime(bin[tls]),bin[l_past],prime(bin[l_past]));
    
    for( int phot=1;phot<Nbin;phot++)
    {   
    sqrt_phot=sqrt(1.*phot);        
    H_fb_l.set(bin[tls](8),prime(bin[tls](4)),bin[r_past](phot),prime(bin[r_past](phot+1)), gam_l_fb*sqrt_phot);
    H_fb_l.set(bin[tls](7),prime(bin[tls](3)),bin[r_past](phot),prime(bin[r_past](phot+1)), gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](6),prime(bin[tls](2)),bin[r_past](phot),prime(bin[r_past](phot+1)), gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](5),prime(bin[tls](1)),bin[r_past](phot),prime(bin[r_past](phot+1)), gam_l_fb*sqrt_phot);        
    // the adjoint left atom interaction hamiltonian
    H_fb_l.set(bin[tls](4),prime(bin[tls](8)),bin[r_past](phot+1),prime(bin[r_past](phot)), c_gam_l_fb*sqrt_phot);
    H_fb_l.set(bin[tls](3),prime(bin[tls](7)),bin[r_past](phot+1),prime(bin[r_past](phot)), c_gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](2),prime(bin[tls](6)),bin[r_past](phot+1),prime(bin[r_past](phot)), c_gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](1),prime(bin[tls](5)),bin[r_past](phot+1),prime(bin[r_past](phot)), c_gam_l_fb*sqrt_phot);                
    // time non-local field ... 
    H_fb_r.set(bin[tls](8),prime(bin[tls](7)),bin[l_past](phot),prime(bin[l_past](phot+1)), gam_r_fb*sqrt_phot);
    H_fb_r.set(bin[tls](6),prime(bin[tls](5)),bin[l_past](phot),prime(bin[l_past](phot+1)), gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](4),prime(bin[tls](3)),bin[l_past](phot),prime(bin[l_past](phot+1)), gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](2),prime(bin[tls](1)),bin[l_past](phot),prime(bin[l_past](phot+1)), gam_r_fb*sqrt_phot);        
    // the adjoint left atom interaction hamiltonian
    H_fb_r.set(bin[tls](7),prime(bin[tls](8)),bin[l_past](phot+1),prime(bin[l_past](phot)), c_gam_r_fb*sqrt_phot);
    H_fb_r.set(bin[tls](5),prime(bin[tls](6)),bin[l_past](phot+1),prime(bin[l_past](phot)), c_gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](3),prime(bin[tls](4)),bin[l_past](phot+1),prime(bin[l_past](phot)), c_gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](1),prime(bin[tls](2)),bin[l_past](phot+1),prime(bin[l_past](phot)), c_gam_r_fb*sqrt_phot);        
    }

    auto H_fb  =   H_fb_l * delta(bin[l_past],prime(bin[l_past])) 
                 + H_fb_r * delta(bin[r_past],prime(bin[r_past]))
               ;
    H_fb  = delta(bin[r_now]   ,prime(bin[r_now]))    * H_fb * delta(bin[l_now]   ,prime(bin[l_now]))      ;
    H_fb  = delta(bin[r_middle],prime(bin[r_middle])) * H_fb * delta(bin[l_middle],prime(bin[l_middle]))   ;
    
    H_dis = delta(bin[l_past],prime(bin[l_past])) * H_dis * delta(bin[r_past],prime(bin[r_past]));
           
    auto H_int = H_dis + H_fb;  //PrintData(H_int);           
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

    auto temp= ITensor(bin[tls],prime(bin[tls]));
    temp.set(bin[tls](1),prime(bin[tls](1)),1.);
    temp.set(bin[tls](2),prime(bin[tls](2)),1.);
    temp.set(bin[tls](3),prime(bin[tls](3)),1.);
    temp.set(bin[tls](4),prime(bin[tls](4)),1.);
    temp.set(bin[tls](5),prime(bin[tls](5)),1.);
    temp.set(bin[tls](6),prime(bin[tls](6)),1.);
    temp.set(bin[tls](7),prime(bin[tls](7)),1.);
    temp.set(bin[tls](8),prime(bin[tls](8)),1.);

    //ITensor temp;
    temp = delta(bin[l_now],prime(bin[l_now]))*temp*delta(bin[l_past],prime(bin[l_past]));
    temp = delta(bin[r_now],prime(bin[r_now]))*temp*delta(bin[r_past],prime(bin[r_past]));
    temp = delta(bin[r_middle],prime(bin[r_middle]))*temp*delta(bin[l_middle],prime(bin[l_middle]));
    // now the taylor expansion of the U_evo

    U_evo =   temp + H_int_1 + H_int_2 + H_int_3 
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
    MPSTensor.set(bin[tls](1),binlink[Nfb](1),binlink[Nfb+1](1),Init[1]); // 000 
    MPSTensor.set(bin[tls](2),binlink[Nfb](1),binlink[Nfb+1](1),Init[2]); // 001
    MPSTensor.set(bin[tls](3),binlink[Nfb](1),binlink[Nfb+1](1),Init[3]); // 010
    MPSTensor.set(bin[tls](4),binlink[Nfb](1),binlink[Nfb+1](1),Init[4]); // 011
    MPSTensor.set(bin[tls](5),binlink[Nfb](1),binlink[Nfb+1](1),Init[5]); // 100
    MPSTensor.set(bin[tls](6),binlink[Nfb](1),binlink[Nfb+1](1),Init[6]); // 101
    MPSTensor.set(bin[tls](7),binlink[Nfb](1),binlink[Nfb+1](1),Init[7]); // 110
    MPSTensor.set(bin[tls](8),binlink[Nfb](1),binlink[Nfb+1](1),Init[8]); // 111
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
    SWAP = psi(k)*psi(k-1);
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


// #################################################################
// #################################################################

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

// #################################################################
// #################################################################

int main(int argc, char* argv[])
{ 
if(argc != 2) { printfln("Usage: %s inputfile",argv[0]); return 0; } 
auto input = InputGroup(argv[1],"input");
// size of timestep
Real dt = input.getReal("time_step");
// maximum number of time steps
int t_end = input.getInt("time_end");
Real Gamma_l = input.getReal("Gamma_l",0.); //Gamma_l = Gamma_l*sqrt(dt);
Real Gamma_m = input.getReal("Gamma_m",0.); //Gamma_m = Gamma_m*sqrt(dt);
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
Init[1] = input.getReal("init_000"); Init[2] = input.getReal("init_001"); 
Init[3] = input.getReal("init_010"); Init[4] = input.getReal("init_011"); 
Init[5] = input.getReal("init_100"); Init[6] = input.getReal("init_101"); 
Init[7] = input.getReal("init_110"); Init[8] = input.getReal("init_111");
    
Real phi_l = input.getReal("phi_l");      
Real phi_m = input.getReal("phi_m");    
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
double pop3=0.;
double prop000 = 
// 000=1, 001=2 ,010=3,011=4,100=5,101=6,110=7,111=8
pop1 = Init[2]*Init[2]+Init[4]*Init[4]+Init[6]*Init[6]+Init[8]*Init[8]; // ATTENTION Init_xyz is real
pop2 = Init[3]*Init[3]+Init[4]*Init[4]+Init[7]*Init[7]+Init[8]*Init[8]; // ATTENTION Init_xyz is real
pop3 = Init[5]*Init[5]+Init[7]*Init[7]+Init[6]*Init[6]+Init[8]*Init[8]; // ATTENTION Init_xyz is real

FILE *file;
FILE *f_prob;

char FILE_NAME[2048+1024];
snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_N3_OCCUPATIONS_FB_at_Step=%i_of%i_sum_pop_%.2f_phase_l_%.2f_phase_m_%.2f_phase_r_%.2f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,fb,t_end,pop3+pop2+pop1,phi_l,phi_m,phi_r,dt,cutoff*1000000.);
file = fopen(FILE_NAME,"w");

snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_N3_PROBABILITY_FB_at_Step=%i_of%i_sum_pop_%.2f_phase_l_%.2f_phase_m_%.2f_phase_r_%.2f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,fb,t_end,pop3+pop2+pop1,phi_l,phi_m,phi_r,dt,cutoff*1000000.);
f_prob = fopen(FILE_NAME,"w");

// this will write the complete config file as a comment in front of the output file
FILE *f_input;
f_input = fopen("parameters.cfg","r");
char line_length [256];
char line_in [256];
while ( fgets(line_length, sizeof line_length, f_input) != NULL )
{ 
  fputs (line_length, stdout); strcpy (line_in, line_length); // getline 
  fprintf(file,"##"); fputs(line_in,file); // write so that grace cannot read it
  fprintf(f_prob,"##"); fputs(line_in,f_prob); // write so that grace cannot read it    
}  
fclose(f_input);
// ----------------------------------------------------------------------------------
// --------------------- SETUP THE MPS ----------------------------------------------
// ----------------------------------------------------------------------------------
// create the physical indices
    
// array starts with zero ... so ttotal + 1 ... ttotal time bins, 1 system bin
auto bin = std::vector<Index>((int)ttotal+1); 
for(int j = 1; j <= ttotal; ++j){ if (j%2==1) bin.at(j) = Index(Nbin,"left");
                                      else        bin.at(j) = Index(Nbin,"right");         }
int tls=0; bin.at(tls) = Index(8,"sys"); 
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
int l_middle = sys_at-Nfb/2   ;
int r_middle = sys_at-Nfb/2 +1;    

ITensor U_evo, temp; 
MPO_SETUP(U_evo,bin,tls,l_past,l_middle,l_now,Nbin,Gamma_l,Gamma_m,Gamma_r,phi_l,phi_m,phi_r,dt);
//PrintData(U_evo);

// ----------------------------------------------------------------------------------
// --------------------- END U_EVO SETUP --------------------------------------------
// ----------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------
// --------------------- START EVOLUTION  -------------------------------------------
// ----------------------------------------------------------------------------------
double mps_norm=1.; mps_norm = TLS_norm(psi(sys_at));
double l_out = 0.;
double r_out= 0.;
double sum_l=0.;
double sum_r=0.;

Real Prob[10]; 
Prob[9]=0.;
for(int p=1;p<=8;p++) {Prob[p]=state_prop(psi(sys_at),p); Prob[9] += Prob[p];}
               
ITensor U,S,V,W,SWAP;
Index iFBl,iFBr,iMl,iMr,iCBl,iCBr; // index fuer feedback, system, und current bin 
Index iLinkp; 
Index b; 
ITensor BdgB;

for(int m=0;m<t_end;m++)
{   
    l_past= sys_at-Nfb    ;
    r_past= sys_at-Nfb+1  ;
    l_now = sys_at        ; 
    r_now = sys_at+1;
    l_middle = sys_at-Nfb/2    ;
    r_middle = sys_at-Nfb/2 +1 ;
    // --- Status and File output ----------------------------------------------------------------------------------------
    sum_l +=l_out; sum_r += r_out;
    printf("Step %i of %i: norm=%.10f  -- ",m,t_end,mps_norm);
    printf("pop_l=%.10f -- pop_m=%.10f -- pop_r=%.10f --",pop3,pop2,pop1); 
    printf("l_out=%.10f -- r_out=%.10f -- Sum_Out=%.2f ",l_out,r_out,sum_l+sum_r); 
    printf("\n");
    fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \t %.10f \n",m*dt,pop3,pop2,pop1,mps_norm);
    fflush(file); 
    fprintf(f_prob,"%.10f \t",m*dt); 
    fprintf(f_prob,"%.10f \t",Prob[1]);
    fprintf(f_prob,"%.10f \t",Prob[2]);fprintf(f_prob,"%.10f \t",Prob[3]);fprintf(f_prob,"%.10f \t",Prob[5]); 
    fprintf(f_prob,"%.10f \t",Prob[4]);fprintf(f_prob,"%.10f \t",Prob[6]);fprintf(f_prob,"%.10f \t",Prob[7]); 
    fprintf(f_prob,"%.10f \t",Prob[8]);
    fprintf(f_prob,"%.10f \t %.10f \t",l_out,r_out);
    fprintf(f_prob,"\n");
    fflush(f_prob);
    // --- End - Status and File output -----------------------------------------------------------------------------------
    //PrintMPS(psi,1,ttotal+1);
    // -------------- PREPARE MPS TO APPLY MPO ------------------------------------------
    if(fb>2)
    { 
      SWAP_RIGHT_FORWARD(psi,bin,r_middle,sys_at-1,cutoff); // move middle bins next to system bin
      SWAP_RIGHT_FORWARD(psi,bin,r_past,sys_at-1,cutoff);   // move beyond middle bins to absorb orthoCenter
    // --------------------------------------------------------
    // --------------------------------------------------------
    // ----------- APPLY MPO START ----------------------------
    iMl = findIndex(psi.A(sys_at-4),"left"); 
    iMr = findIndex(psi.A(sys_at-3),"right");      
    iFBl = findIndex(psi.A(sys_at-2),"left"); 
    iFBr = findIndex(psi.A(sys_at-1),"right");      
    iCBl = findIndex(psi.A(sys_at+1),"left");
    iCBr = findIndex(psi.A(sys_at+2),"right");    
    iLinkp = commonIndex(psi.A(sys_at-4),psi.A(sys_at-5));
    
    temp = noPrime(psi(sys_at-4)*psi(sys_at-3)*psi(sys_at-2)*psi(sys_at-1)*U_evo*psi(sys_at)*psi(sys_at+1)*psi(sys_at+2));  
    //Print(temp);
    
    U=ITensor(iLinkp,iMl,iMr,iFBl,iFBr,iCBl,iCBr); 
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    // tls in the next mps slot 
    psi.setA(sys_at+2,V); 
    Prob[9]=0;
    for(int p=1;p<=8;p++) {Prob[p]=state_prop(V*S,p); Prob[9] += Prob[p]; }
    pop1 =  Prob[2]+Prob[4]+Prob[6]+Prob[8];
    pop2 =  Prob[3]+Prob[4]+Prob[7]+Prob[8];
    pop3 =  Prob[5]+Prob[6]+Prob[7]+Prob[8];  
    mps_norm = eltC( dag(V*S)*V*S).real(); 
    // now factorize step by step the reservoir tensor moving orthoCenter
    temp=U*S; U=ITensor(iLinkp,iMl,iMr,iFBl,iFBr,iCBl); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at+1,V); // iCBr 
    temp=U*S; U=ITensor(iLinkp,iMl,iMr,iFBl,iFBr); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at,V); // iCBl 
    temp=U*S; U=ITensor(iLinkp,iMl,iMr,iFBl); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at-1,V); // iFBr
    // ---------------------------------------
    b=findIndex(V,"right");
    BdgB = ITensor(b,prime(b)); 
    BdgB.set(b(2),prime(b(2)),1.); // ATTENTION bath bin maximal 1 photon!!! Nbin 2!!! 
    r_out = eltC( dag(V*S)*noPrime(BdgB*V*S) ).real();
    // ---------------------------------------
    temp=U*S; U=ITensor(iLinkp,iMl,iMr); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at-2,V); // iFBl 
    // ---------------------------------------
    b=findIndex(V,"left");
    BdgB = ITensor(b,prime(b)); 
    BdgB.set(b(2),prime(b(2)),1.); // ATTENTION bath bin maximal 1 photon!!! Nbin 2!!! 
    l_out = eltC( dag(V*S) * noPrime(BdgB*V*S) ).real();
    // ---------------------------------------
    temp=U*S; U=ITensor(iLinkp,iMl); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at-3,V); // iMr 
    psi.setA(sys_at-4,U*S); //iMl    
    // ----------- SWAP BACK START -----------------------------
    SWAP_LEFT_BACKWARD(psi,bin,sys_at-4,l_middle-2,cutoff); 
    SWAP_LEFT_BACKWARD(psi,bin,sys_at-2,l_past,cutoff); 
    // ----------- SWAP BACK END ------------------------------
    // --------------------------------------------------------
    }
    else
    {
    iFBl = findIndex(psi.A(sys_at-4),"left"); 
    iFBr = findIndex(psi.A(sys_at-3),"right");      
    iMl = findIndex(psi.A(sys_at-2),"left"); 
    iMr = findIndex(psi.A(sys_at-1),"right");      
    iCBl = findIndex(psi.A(sys_at+1),"left");
    iCBr = findIndex(psi.A(sys_at+2),"right");    
    iLinkp = commonIndex(psi.A(sys_at-4),psi.A(sys_at-5));
    
    temp = noPrime(psi(sys_at-4)*psi(sys_at-3)*psi(sys_at-2)*psi(sys_at-1)*U_evo*psi(sys_at)*psi(sys_at+1)*psi(sys_at+2));  
    
    U=ITensor(iLinkp,iMl,iMr,iFBl,iFBr,iCBl,iCBr);     
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    // tls in the next mps slot 
    psi.setA(sys_at+2,V); 
    Prob[9]=0;
    for(int p=1;p<=8;p++) {Prob[p]=state_prop(V*S,p); Prob[9] += Prob[p]; }
    pop1 =  Prob[2]+Prob[4]+Prob[6]+Prob[8];
    pop2 =  Prob[3]+Prob[4]+Prob[7]+Prob[8];
    pop3 =  Prob[5]+Prob[6]+Prob[7]+Prob[8];  
    mps_norm = eltC( dag(V*S)*V*S).real(); 
    // now factorize step by step the reservoir tensor moving orthoCenter
    temp=U*S; U=ITensor(iLinkp,iFBl,iFBr,iMl,iMr,iCBl); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at+1,V); // iCBr 
    temp=U*S; U=ITensor(iLinkp,iFBl,iFBr,iMl,iMr); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at,V); // iCBl 
    temp=U*S; U=ITensor(iLinkp,iFBl,iFBr,iMl); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at-1,V); // iMr
    temp=U*S; U=ITensor(iLinkp,iFBl,iFBr); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at-2,V); // iMl 
    temp=U*S; U=ITensor(iLinkp,iFBl); svd(temp,U,S,V,{"Cutoff=",cutoff}); 
    // ---------------------------------------
    b=findIndex(V,"right");
    BdgB = ITensor(b,prime(b)); 
    BdgB.set(b(2),prime(b(2)),1.); // ATTENTION bath bin maximal 1 photon!!! Nbin 2!!! 
    r_out = eltC( dag(V*S)*noPrime(BdgB*V*S) ).real();
    // ---------------------------------------
    psi.setA(sys_at-3,V*S); // iFBr 
    // ---------------------------------------
    b=findIndex(U,"left");
    BdgB = ITensor(b,prime(b)); 
    BdgB.set(b(2),prime(b(2)),1.); // ATTENTION bath bin maximal 1 photon!!! Nbin 2!!! 
    l_out = eltC( dag(U*S) * noPrime(BdgB*U*S) ).real();
    // ---------------------------------------
    psi.setA(sys_at-4,U); //iFBl      
    
    temp = psi(sys_at-3)*psi(sys_at-2);
    iLinkp=commonIndex(psi(sys_at-3),psi(sys_at-4));  
    U = ITensor(iFBr,iLinkp);
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    psi.setA(sys_at-3,U);       
    psi.setA(sys_at-2,V*S);      
    }
    // Update MPO, only if a further loop is possible
    if (m<t_end-1)
    {       
    U_evo = delta(bin[l_now],bin[l_now+2])*U_evo*delta(prime(bin[l_now]),prime(bin[l_now+2])) ;
    U_evo = delta(bin[r_now],bin[r_now+2])*U_evo*delta(prime(bin[r_now]),prime(bin[r_now+2])) ;
    U_evo = delta(bin[l_middle],bin[l_middle+2])*U_evo*delta(prime(bin[l_middle]),prime(bin[l_middle+2])) ;
    U_evo = delta(bin[r_middle],bin[r_middle+2])*U_evo*delta(prime(bin[r_middle]),prime(bin[r_middle+2])) ;
    U_evo = delta(bin[l_past],bin[l_past+2])*U_evo*delta(prime(bin[l_past]),prime(bin[l_past+2])) ;
    U_evo = delta(bin[r_past],bin[r_past+2])*U_evo*delta(prime(bin[r_past]),prime(bin[r_past+2])) ;
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
