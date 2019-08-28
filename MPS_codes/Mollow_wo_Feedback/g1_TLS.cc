#include "itensor/all.h"

using namespace itensor;

double population(ITensor A, int Nbin)
{
    Index s=findIndex(A,"sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(2),prime(s(2)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).real();
}

double coherence(ITensor A, int Nbin)
{
    Index s=findIndex(A,"sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(1),prime(s(2)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).imag();
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
    //timestep
    Real dt = input.getReal("time_step");
    //end of integration
    Real t_end = input.getReal("time_end");
    //number of steps
    Real ttotal = t_end/dt;
    //coherent pumping strengths
    Real Omega = input.getReal("Omega_TLS",0.);
    //TLS decay rate
    Real Gamma = input.getReal("Gamma",0.);
    //dimension of local Hilbert space of each bin
    int Nbin = input.getInt("Nbin",4); 
    //cutoff of schmidt values
    Real cutoff = input.getReal("svdcutoff");
    //maximal number of schmidtvalues 
    int maxm = input.getInt("maxnumberofSV");
    Real init_EE = input.getReal("init_EE");
    Real init_GG = input.getReal("init_GG");
    Real benchmark_factor = input.getReal("benchmark_factor");
    int SpectrumSteps = input.getInt("SpectrumSteps");
    Real SpectrumIntervall = input.getReal("SpectrumIntervall");
    FILE * file;
    file = fopen ("system_dynamics_output.dat","w");
    fprintf(file,"#time \tS+S- \t norm \n");
    //dimension of system bin
    int d = 2; 
    int gtotal =2*((int)ttotal);

    // ###################################################################################
    // setup physical degrees of freedom bins from 1 to gtotal+1
    auto bin = std::vector<Index>(gtotal+2);
    bin.at(1) = Index(Nbin,"sys");
    for(int j = 2; j <= gtotal+1; ++j){ bin.at(j) = Index(Nbin,"bath"); }
    // setup links in between the physical degrees of freedom
    auto binlink = std::vector<Index>(gtotal+1);
    for(int j = 1; j <= gtotal; ++j) { binlink.at(j) = Index(1,"link");   }
    // ###################################################################################
    // now initialize an MPS 
    MPS psi   = MPS(gtotal+1);
    ITensor MPSTensor = ITensor();
    
    // first bin loaded with system and its initial conditions
    MPSTensor = ITensor(bin[1],binlink[1]);
    MPSTensor.set(bin[1](2),binlink[1](1),init_EE);
    MPSTensor.set(bin[1](1),binlink[1](1),init_GG);
    psi.setA(1,MPSTensor);
    // bins in between 1 and gototal+1 have links left and right and are in the ground state      
    for(int j = 2; j<=gtotal; j++)
    {
      MPSTensor = ITensor(bin[j],binlink[j-1],binlink[j]);
      MPSTensor.set(bin[j](1),binlink[j-1](1),binlink[j](1),1.);
      psi.setA((j),MPSTensor);
    }
    // last bin has links only to the left
    MPSTensor = ITensor(bin[gtotal+1],binlink[gtotal]);
    MPSTensor.set(bin[gtotal+1](1),binlink[gtotal](1),1.);
    psi.setA((gtotal+1),MPSTensor);
    // now I have an MPS, filled with bin[1] bin[2] ... bin[gtotal+1]
    // bin[1] is the system, and bin[1] is changed throughout the algorithm
    // set orthoCenter
    //psi.position(1);
    psi.orthogonalize();
    //printf("%i \n",orthoCenter(psi));

   // ###################################################################################
   // allocate system evo matrix    
    auto H_sys  = ITensor(bin[1],prime(bin[1]));
    //allocate H_int evo matrix - interaction of system with first timebin
    auto H_dis = ITensor(bin[1],prime(bin[1]),bin[2],prime(bin[2]));
        
    //pumped TLS 
    H_sys.set(bin[1](1),prime(bin[1](2)),-Cplx_i*dt*Omega);
    H_sys.set(bin[1](2),prime(bin[1](1)),-Cplx_i*dt*Omega);
    
    //H_int with photon creation/destruction in bin k
    H_dis.set(bin[1](2),prime(bin[1](1)),bin[2](1),prime(bin[2](2)),-Cplx_i*sqrt(dt)*Gamma);
    H_dis.set(bin[1](1),prime(bin[1](2)),bin[2](2),prime(bin[2](1)),-Cplx_i*sqrt(dt)*Gamma);
       
    auto H_int = H_sys*delta(bin[2],prime(bin[2])) + H_dis ;
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
 
    auto temp= ITensor(bin[1],prime(bin[1]));
    temp.set(bin[1](1),prime(bin[1](1)),1.);
    temp.set(bin[1](2),prime(bin[1](2)),1.);
    
    auto U_evo = delta(bin[2],prime(bin[2]))*temp
                 + H_int_1 + H_int_2 + H_int_3 + H_int_4 + H_int_5 
                 + H_int_6 + H_int_7 + H_int_8 + H_int_9 + H_int_10 
               ;
    
    // ###################################################################################

    // ###################################################################################
    ITensor U,S,Vp;
    int j;

    double rhoEE=init_EE;
    double rhoGE=0.;
 
    // variables for the exact solution obtained via Laplace transform
    double tGamma = benchmark_factor*Gamma*Gamma;///(2.*3.141);
    double tOmega = 2.*Omega;
    double bOm =sqrt(1.*tOmega*tOmega-tGamma*tGamma/16.);
    double bGam=3.*tGamma/4.;    
    double exact;
    double eG;
    double sOm;
    double cOm;
    double S0 =init_EE+init_GG;
    double W0 =init_EE-init_GG;
    printf("dt=%.10f \n",dt);
    fprintf(file,"%.10f  \t %.10f \t %.10f \t %.10f \t %.10f \n",0.,rhoEE,init_EE,0.,1.0000000);
    
    for(j = 1; j<gtotal-10;j++)
    {
        if(j>1){ U_evo=U_evo*delta(bin[j],bin[j+1])*delta(prime(bin[j]),prime(bin[j+1]));  }

        // fuer psi ####################################################
        temp=noPrime(U_evo * psi.A(j) * psi.A(j+1));
        if (j==1) U=ITensor(bin[j+1]); 
        else U=ITensor(bin[j+1],commonIndex(psi.A(j),psi.A(j-1)));
        svd(temp,U,S,Vp);
        Vp=Vp*S;       
        psi.setA(j,U);
        psi.setA(j+1,Vp);
        psi.position(j+1);
        
        // exact solution of the Mollow problem without initial polarization
        sOm  = sin(dt*j*bOm);
        cOm  = cos(dt*j*bOm);
        eG   = exp(-dt*j*bGam);
        
        if (Omega>0.0000001)
        {    
        exact  = eG*sOm*( -S0     *tGamma/bOm
                          -W0*0.25*tGamma/bOm
                          +S0     *(bGam/bOm)*tGamma*tGamma/(2.*tOmega*tOmega+tGamma*tGamma) );
        exact += eG*cOm*(  S0                *tGamma*tGamma/(2.*tOmega*tOmega+tGamma*tGamma) 
                         + W0 ) ;
        exact +=          -S0                *tGamma*tGamma/(2.*tOmega*tOmega+tGamma*tGamma);
        }
        else exact = -S0+2.*init_EE*exp(-tGamma*j*dt); 

        if ( gtotal > 1000 ) // if number of time steps too large, just save every 1000th data points
        {
        if    ((j % (gtotal/1000)) == 0) 
        { fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \t %.10f \n",j*dt,population(psi(j+1),Nbin),0.5*exact+0.5*S0,coherence(psi(j+1),Nbin),norm(psi)); }
        } 
        else 
        { 
        fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \t %.10f \n",j*dt,population(psi(j+1),Nbin),0.5*exact+0.5*S0,coherence(psi(j+1),Nbin),norm(psi)); 
        }
        
        if  ( ((j % (gtotal/20)) == 0) || (j>gtotal+20) ) 
        {  printf("%f \t %f---%f \t %f -- j=%i \n",j*dt,population(psi(j+1),Nbin),0.5*exact+0.5*S0,coherence(psi(j+1),Nbin),j);                        }
                   
    }
    
    fclose(file);
    // ################################################################################### 
    file = fopen ("bath_dynamics_output.dat","w");
    // now we have a steady state and can calculate the spectrum <b^\dg(t)b(t-\tau)>
    // setup the bath correlation array
    auto bath_g1 = std::vector<Complex>(SpectrumSteps+1);
    auto bath_g2= std::vector<Complex>(SpectrumSteps+1);
    ITensor Tg1;
    ITensor Tg2;
    
    // ------------ swap ortoCenter to the bath bin -----------------
    temp=psi(j)*psi(j-1); //PrintData(temp);
    int RefBin=gtotal-10;
    U = ITensor(bin[RefBin],commonIndex(psi(j-2),psi(j-1))); //PrintData(U);
    svd(temp,U,S,Tg1);
    psi.setA(j-1,U*S);
    psi.setA(j,Tg1); 
    // --------- end prepare via swap -------------------------------

    // now define the flip operator for RefBin
    Index s=findIndex(psi(j-1),"bath"); //Print(s);
    ITensor Sp = ITensor(s,prime(s));
    Sp.set(s(2),prime(s(1)),1.); //PrintData(Sp);

    // first we need the steady state value of the bath coherence
     Tg1= dag( psi(j-1) )* noPrime(Sp*psi(j-1));
     Complex bc_inf = eltC(Tg1);
     double  b_offset = abs(bc_inf)*abs(bc_inf);
     Tg2= dag( noPrime(Sp*psi(j-1)) ) * noPrime(noPrime(Sp*psi(j-1)));
     double bdgb = eltC(Tg2).real(); 
     printf("bath_coherence=%.10f+I%.10f -- |b_c|^2=%.10f -- <b^+b>=%.10f \n",bc_inf.real(),bc_inf.imag(),b_offset,bdgb);
    
    // the first entry is the self-correlation of the reference bin
     Tg1= dag( noPrime(psi(j-1)*Sp) )*noPrime(Sp*psi(j-1));
     Tg2= dag( psi(j-1) )* noPrime(Sp*noPrime(Sp*psi(j-1)));
     bath_g1[0]=eltC(Tg1);
     bath_g2[0]=eltC(Tg2);
     
     printf("bath_g1[%i]-|bc_inf|^2=%.10f+I%.10f -- bath_g2[0]=%.10f \n",0,bath_g1[0].real()-b_offset,bath_g1[0].imag(),bath_g2[0].real());


     fprintf(file,"\n");
     fprintf(file,"%.10f \t %.10f \t %10f \t %.10f \n",dt*0,bath_g1[0].real()-b_offset,bath_g1[0].imag(),bath_g2[0].real());
    //Index for CorrelatedBin
    Index t=findIndex(psi(j-2),"bath"); //Print(t);
    ITensor Sm = ITensor(t,prime(t));
    Sm.set(t(2),prime(t(1)),1.);
    
    for(int i=j-1;i>j-SpectrumSteps;i--)
    {   
    t=findIndex(psi(i-1),"bath"); //Print(t);
    Sm = ITensor(t,prime(t));
    Sm.set(t(2),prime(t(1)),1.);

    temp = psi(i)*psi(i-1);

    Tg1= dag( noPrime(temp*Sp) )*noPrime(Sm*temp);
    Tg2= dag( noPrime(Sm*noPrime(Sp*temp) )) * noPrime(Sm*noPrime(Sp*temp));
    bath_g1[j-i+1]=eltC(Tg1);
    bath_g2[j-i+1]=eltC(Tg2);
    printf("bath_g1[%i]-|bc_inf|^2=%.10f+I%.10f -- g2=%.10f \n",j-i+1,bath_g1[j-i+1].real()-b_offset,bath_g1[j-i+1].imag(),bath_g2[j-i+1].real());
    fprintf(file,"%.10f \t %.10f \t %10f \t %.10f \n",dt*(j-i+1),bath_g1[j-i+1].real()-b_offset,bath_g1[j-i+1].imag(),bath_g2[j-i+1].real());

    U = ITensor(bin[RefBin],commonIndex(psi(i-1),psi(i-2))); //PrintData(U);
    svd(temp,U,S,Tg1);
    psi.setA(i-1,U*S);
    psi.setA(i,Tg1); 
    }

    fclose(file);    
// now we have the bath correlation, saved into the array bath_g1, given these correlations, the spectrum
// can be computed
   
   file = fopen ("spectrum_g2.dat","w");

   auto spectrum = std::vector<Complex>(SpectrumSteps+1);
   double dw = SpectrumIntervall/SpectrumSteps; 
   double om;
   
   for(int i=0;i<=SpectrumSteps;i++)
   {
     om = -0.5*(SpectrumIntervall)+i*dw;
     spectrum[i] = 0.;
     for(int b=0;b<=SpectrumSteps;b++) spectrum[i] += exp(-Cplx_i*om*b*dt)*bath_g1[b];  
     fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \n",om,spectrum[i].real(),spectrum[i].imag(),bath_g2[i].real()/(bdgb*bdgb));
   }

   fclose(file);
   
return 0;        
}
