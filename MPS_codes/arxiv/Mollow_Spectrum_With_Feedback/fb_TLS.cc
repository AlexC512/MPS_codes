#include "itensor/all.h"

using namespace itensor;

double TLS_occupation(ITensor A)
{
    //search for Index in Tensor (should be system)
    Index s=findIndex(A, "sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(2),prime(s(2)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).real(); 
}

double TLS_norm(ITensor A)
{
    return eltC( dag(A)*A).real(); 
}

/*
void SWAP_FORWARD(MPS psi,std::vector<Index> bin, int i, int Nfb, double cutoff)
{
    ITensor SWAP,U,S,V;
//    PrintData(psi(1));PrintData(psi(2));
    for(int k=1;k<Nfb-1;k++) // swap feedback bin next to tls bin
    {    
            SWAP = psi.A(i-Nfb+k)*psi.A(i-Nfb+k+1); 
            if (i-Nfb+k==1) U=ITensor(bin[i-Nfb+k+1]); // keine Linksverlinkung fuer psi.A(1)
            else U=ITensor(bin[i-Nfb+k+1],commonIndex(psi.A(i-Nfb+k),psi.A(i-Nfb+k-1)));  
            svd(SWAP,U,S,V,{"Cutoff=",cutoff});
            psi.setA(i-Nfb+k,U); 
            psi.setA(i-Nfb+k+1,S*V);
    }

//    PrintData(psi(1));
 return;
} */
    
int main(int argc, char* argv[])
    { 
        
time_t curtime; struct tm *loctime; curtime = time(NULL); loctime = localtime (&curtime); 
int hour = loctime -> tm_hour;int minute = loctime -> tm_min;int second = loctime -> tm_sec;
int year = 1900+loctime -> tm_year;int month = loctime -> tm_mon;int day  = loctime -> tm_mday;
printf("Date: %d.%d.%.d -- Time: %d:%d:%d  \n",day,month+1,year,hour,minute,second);  
        
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
    int Nfb = input.getInt("feedback_time");
    Real init_22 = input.getReal("init_excited_state"); 
    Real init_11 = input.getReal("init_ground_state"); 
    Real phi = input.getReal("phi"); 
    Real fb_on = input.getReal("fb_on");
    int SpectrumSteps = input.getInt("SpectrumSteps");
    Real SpectrumIntervall = input.getReal("SpectrumIntervall");    
    
    printf("TLS initial state: (22)=%.2f -- (11)=%.2f \n",init_22,init_11);
    if (Nfb == 0) Nfb = ttotal;
    //dimension of system bin
    int d = 2; 

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // --------------------- SETUP THE MPS ----------------------------------------------
    // ----------------------------------------------------------------------------------
    // create the physical indices
    auto bin = std::vector<Index>((int)ttotal+3);
    for(int j = 1; j <= ttotal+1; ++j){ bin.at(j) = Index(Nbin,"bath"); }
    int tls=ttotal+2;
    bin.at(tls) = Index(Nbin,"sys"); // bin(1) sys throughout the code

    auto binlink = std::vector<Index>((int)ttotal+1);
    for(int j = 1; j <= ttotal; ++j) { binlink.at(j) = Index(1,"Link");   }
    MPS psi=MPS(ttotal+1);
    ITensor MPSTensor = ITensor();

    // set up the mps 
    // first MPS slot is filled with a bin which has no link to the left
     MPSTensor = ITensor(bin[1],binlink[1]);
     MPSTensor.set(bin[1](1),binlink[1](1),1.); // Level besetzt im TLS
     psi.setA((1),MPSTensor);

    for(int j = 2; j<Nfb; j++)
    {
     MPSTensor = ITensor(bin[j],binlink[j-1],binlink[j]);
     MPSTensor.set(bin[j](1),binlink[j-1](1),binlink[j](1),1.); // Level besetzt im TLS
     psi.setA((j),MPSTensor);
    }   
    // at MPS slot Nfb the tls tensor is put
    MPSTensor = ITensor(bin[tls],binlink[Nfb-1],binlink[Nfb]);
    MPSTensor.set(bin[tls](2),binlink[Nfb-1](1),binlink[Nfb](1),init_22); // Level besetzt im TLS
    MPSTensor.set(bin[tls](1),binlink[Nfb-1](1),binlink[Nfb](1),init_11); // Level besetzt im TLS
    psi.setA((Nfb),MPSTensor);
    
    for(int j = Nfb; j<=ttotal-1; j++)
    {
     MPSTensor = ITensor(bin[j],binlink[j],binlink[j+1]);
     MPSTensor.set(bin[j](1),binlink[j](1),binlink[j+1](1),1.); // Level besetzt im TLS
     psi.setA((j+1),MPSTensor);
    }
    // the last slot of MPS is filled with a tensor with only links to the left
    MPSTensor = ITensor(bin[ttotal+1],binlink[ttotal]);
    MPSTensor.set(bin[ttotal+1](1),binlink[ttotal](1),1.); // Level besetzt im TLS
    psi.setA((ttotal+1),MPSTensor);
    psi.position(1);
    printf("OrthoCenter=%i \n",orthoCenter(psi));

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ---------------------END MPS SETUP -----------------------------------------------
    // ----------------------------------------------------------------------------------

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ---------------------U_EVO SETUP -------------------------------------------------
    // ----------------------------------------------------------------------------------
    auto U_sysfb =ITensor(bin[tls],prime(bin[tls]));
    // interaction with futurebin
    auto U_disfb = ITensor(bin[tls],prime(bin[tls]),bin[Nfb],prime(bin[Nfb]));
    // interaction with the pastbin
    auto U_fb    = ITensor(bin[tls],prime(bin[tls]),bin[1],prime(bin[1]));
        
    //pumped TLS 
    U_sysfb.set(bin[tls](1),prime(bin[tls](2)),-Cplx_i*dt*Omega);
    U_sysfb.set(bin[tls](2),prime(bin[tls](1)),-Cplx_i*dt*Omega);
    
    //H_int with photon creation/destruction in bin k
   for(int phot=1;phot<Nbin;phot++)
   {            
    U_disfb.set(bin[tls](2),prime(bin[tls](1)),bin[Nfb](phot  ),prime(bin[Nfb](phot+1)),( 1.)*sqrt(dt)*Gamma*sqrt(phot));
    U_disfb.set(bin[tls](1),prime(bin[tls](2)),bin[Nfb](phot+1),prime(bin[Nfb](phot  )),(-1.)*sqrt(dt)*Gamma*sqrt(phot));        
    U_fb.set(bin[tls](2),prime(bin[tls](1)),bin[1](phot  ),prime(bin[1](phot+1)),(-1.)*sqrt(dt)*Gamma*exp(-Cplx_i*phi)*sqrt(phot));
    U_fb.set(bin[tls](1),prime(bin[tls](2)),bin[1](phot+1),prime(bin[1](phot  )),( 1.)*sqrt(dt)*Gamma*exp(Cplx_i*phi)*sqrt(phot));
   }//phot

    auto U_int = U_disfb*delta(bin[1],prime(bin[1])) 
                +U_fb   *delta(bin[Nfb],prime(bin[Nfb]))    ;
    
    ITensor temp;            
    auto tempfb = ITensor(bin[tls],prime(bin[tls]));
    tempfb.set(bin[tls](1),prime(bin[tls](1)),1.);
    tempfb.set(bin[tls](2),prime(bin[tls](2)),1.);

    auto U_fb_1  = U_sysfb*delta(bin[Nfb],prime(bin[Nfb]))*delta(bin[1],prime(bin[1])) +U_int; 
    auto U_fb_2  = (1./2.)  * mapPrime(U_fb_1*prime(U_fb_1),2,1);
    auto U_fb_3  = (1./3.)  * mapPrime(U_fb_1*prime(U_fb_2),2,1);
    auto U_fb_4  = (1./4.)  * mapPrime(U_fb_1*prime(U_fb_3),2,1);
    auto U_fb_5  = (1./5.)  * mapPrime(U_fb_1*prime(U_fb_4),2,1);
    auto U_fb_6  = (1./6.)  * mapPrime(U_fb_1*prime(U_fb_5),2,1);
    auto U_fb_7  = (1./7.)  * mapPrime(U_fb_1*prime(U_fb_6),2,1);
    auto U_fb_8  = (1./8.)  * mapPrime(U_fb_1*prime(U_fb_7),2,1);
    auto U_fb_9  = (1./9.)  * mapPrime(U_fb_1*prime(U_fb_8),2,1);
    auto U_fb_10 = (1./10.) * mapPrime(U_fb_1*prime(U_fb_9),2,1);
    
    auto U_evofb =   delta(bin[Nfb],prime(bin[Nfb]))*tempfb*delta(bin[1],prime(bin[1])) 
                   + U_fb_1 + U_fb_2 + U_fb_3 + U_fb_4 + U_fb_5 + U_fb_6 + U_fb_7  + U_fb_8 + U_fb_9 + U_fb_10 ;
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // --------------------- END U_EVO SETUP --------------------------------------------
    // ----------------------------------------------------------------------------------

    
    ITensor U,S,V,W,SWAP;
    int j=1;
    
    double cv_norm,cv_pop;
    
// das bild im MPS sieht aus [pastbin][feedbackbin][systembin][currentbin]
Index iFB,iSB,iCB,iPB; // index fuer feedback, system, und current bin 
Index iFBPB; // index fuer feedback bin und past bin
// zum zeit berechnen
    FILE * file;
    char FILE_NAME[2048+1024];
    snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_TLS_w_FEEDBACK_at_Step=%i_of%i_OmL_%.3f_Gamma%.3f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,Nfb,(int)ttotal,Omega,Gamma,dt,cutoff*1000000.);
    file = fopen(FILE_NAME,"w");
    fprintf(file,"#time \tS+S- \t norm \n");
    fprintf(file,"# %02d_%02d_%02d_%02d_%02d_%02d \n",year,month+1,day,hour,minute,second);
    fprintf(file,"%f \t %e \t %e \n",0.,init_22,0.0000);
    

int ehour = loctime -> tm_hour; int eminute = loctime -> tm_min; int esecond = loctime -> tm_sec;
int eyear = 1900+loctime -> tm_year; int emonth = loctime -> tm_mon; int eday  = loctime -> tm_mday;
double time_before,time_now; time_before =esecond + 60*eminute+60*60*ehour;

for(int i=Nfb;i<ttotal;i++)            
{   
    for(int k=1;k<Nfb-1;k++) // swap feedback bin next to tls bin
    {
            SWAP = psi.A(i-Nfb+k)*psi.A(i-Nfb+k+1); 
            if (i-Nfb+k==1) U=ITensor(bin[i-Nfb+k+1]); // keine Linksverlinkung fuer psi.A(1)
            else U=ITensor(bin[i-Nfb+k+1],commonIndex(psi.A(i-Nfb+k),psi.A(i-Nfb+k-1)));  
            svd(SWAP,U,S,V,{"Cutoff=",cutoff});
            psi.setA(i-Nfb+k,U); 
            psi.setA(i-Nfb+k+1,S*V);
    }
      // define the links for feedback bin, sys bin and current bin
     iFB = findIndex(psi.A(i-1),"bath"); iFBPB = commonIndex(psi.A(i-1),psi.A(i-2));
     iSB = findIndex(psi.A(i),"sys");
     iCB = findIndex(psi.A(i+1),"bath");
     // now apply the mpo
     temp=noPrime(U_evofb * psi.A(i-1) * psi.A(i) * psi.A(i+1)); 
     U=ITensor(iCB,iFB,iFBPB); 
     svd(temp,U,S,V,{"Cutoff=",cutoff});
     // tls in the next mps slot 
     psi.setA(i+1,V); 
     cv_norm = TLS_norm(V*S); 
     cv_pop  = TLS_occupation(V*S); // no access to population directly ... because of link indices ... tensor of 3rd order
     fprintf(file,"%.10f \t %.10f \t %10.f \n",(i-Nfb+1)*dt,cv_pop,cv_norm);
     // now factorize feedback bin and current bin
     W=U*S; 
     U=ITensor(iFB,iFBPB); 
     svd(W,U,S,V,{"Cutoff=",cutoff});
     psi.setA(i-1,U*S); 
     psi.setA(i,V);
     // bring feedback bin to its original position in the mps -- swap
     for(int k=1;k<Nfb-1;k++)
     {      
            SWAP = psi.A(i-1-k)*psi.A(i-k);
            if ( (i-k-1)==1 ) U=ITensor(bin[i-Nfb+1]); // kein zu vererbender Link
            else U=ITensor(bin[i-Nfb+1],commonIndex(psi.A(i-k-2),psi.A(i-k-1)));
            svd(SWAP,U,S,V,{"Cutoff=",cutoff}); 
            // -- lasse Orthocenter beim vorletzten Bin, der im naechsten Schritt wieder zum Systembin geswappt wird
            if (k==Nfb-2) { psi.setA(i-1-k,U);       psi.setA(i-k,V*S); }
            else          { psi.setA(i-1-k,U*S);       psi.setA(i-k,V); }
     }
     if (i<ttotal-1)
     {    
     // define new evolution operator by replacing the old feedback bin with the next one and a new current bin ... tls bin stays
     // however only if there is a next time step
     U_evofb=U_evofb*delta(bin[i],bin[i+1])*delta(prime(bin[i]),prime(bin[i+1]))*delta(bin[i-Nfb+1],bin[i-Nfb+2])*delta(prime(bin[i-Nfb+1]),prime(bin[i-Nfb+2]));
     }
     // show the calculation status every percent
     if ( (i % ((int)ttotal/100)) == 0) 
       { 
         curtime = time(NULL); loctime = localtime (&curtime); 
         ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
         time_now = esecond + 60*eminute+60*60*ehour;  
         printf("Feedback time %.3f of %.3f -- pop=%.10f -- norm=%.10f -- Time elapsed=%.2f sec. \n",i*dt,t_end*1.,cv_pop,cv_norm,time_now-time_before);
         time_before = time_now;
       }     
     fflush(file); // This will flush any pending fprintf output
}//end timeloop

// give information how long it took
curtime = time(NULL); loctime = localtime (&curtime); 
ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
eyear = 1900+loctime -> tm_year; emonth = loctime -> tm_mon; eday  = loctime -> tm_mday;
printf("End of calculation: %d.%d.%.d -- Time: %d:%d:%d  \n",eday,emonth+1,eyear,ehour,eminute,esecond);  
fprintf(file,"# end of calculation %02d_%02d_%02d_%02d_%02d_%02d \n",eyear,emonth+1,eday,ehour,eminute,esecond);
printf("Time elapsed: %d.%d.%.d -- Time: %dh:%dmin:%dsec \n",eday-day,emonth-month,eyear-year,ehour-hour,eminute-minute,esecond-second);  
fprintf(file,"# time difference %02d_%02d_%02d_%02d_%02d_%02d \n",eyear-year,emonth-month,eday-day,ehour-hour,eminute-minute,esecond-second);
// close the output file pointer    
fclose(file);

// observables .... 


//PrintData(psi(ttotal-Nfb-1));
//PrintData(psi(ttotal-Nfb));
//PrintData(psi(ttotal-Nfb+1)); // systembin

temp = dag(psi(ttotal-Nfb+1))*psi(ttotal-Nfb+1); // this is the next feedback bin with OrthoCenter 
PrintData(temp); //swap orthoCenter to the previous one

temp = psi(ttotal-Nfb+1)*psi(ttotal-Nfb);
U = ITensor(bin[ttotal-Nfb+1],commonIndex(psi(ttotal-Nfb+1),psi(ttotal-Nfb+2)));
svd(temp,U,S,V);
psi.setA(ttotal-Nfb+1,U);
psi.setA(ttotal-Nfb,V*S);

temp = dag(psi(ttotal-Nfb))*psi(ttotal-Nfb); // now the first non-interacting past bin is normalized
PrintData(temp);

int RefBin=ttotal-Nfb;

snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_TLS_BathCorrelations_w_FEEDBACK_at_Step=%i_of%i_OmL_%.3f_Gamma%.3f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,Nfb,(int)ttotal,Omega,Gamma,dt,cutoff*1000000.);
    file = fopen(FILE_NAME,"w");
    // now we have a steady state and can calculate the spectrum <b^\dg(t)b(t-\tau)>
    // setup the bath correlation array
    auto bath_g1 = std::vector<Complex>(SpectrumSteps+1);
    auto bath_g2= std::vector<Complex>(SpectrumSteps+1);
    ITensor Tg1;
    ITensor Tg2;
    
    // now define the flip operator for RefBin
    Index s=findIndex(psi(RefBin),"bath"); //Print(s);
    ITensor Sp = ITensor(s,prime(s));
    Sp.set(s(2),prime(s(1)),1.); //PrintData(Sp);
    
    // the first entry is the self-correlation of the reference bin
     Tg1= dag( noPrime(psi(RefBin)*Sp) )*noPrime(Sp*psi(RefBin));
     Tg2= dag( psi(RefBin) )* noPrime(Sp*noPrime(Sp*psi(RefBin)));
     bath_g1[0]=eltC(Tg1);
     bath_g2[0]=eltC(Tg2);
     
     printf("bath_g1[%i]-|bc_inf|^2=%.10f+I%.10f -- bath_g2[0]=%.10f \n",0,bath_g1[0].real(),bath_g1[0].imag(),bath_g2[0].real());

     fprintf(file,"%.10f \t %.10f \t %10f \t %.10f \n",dt*0,bath_g1[0].real(),bath_g1[0].imag(),bath_g2[0].real());
    //Index for CorrelatedBin
    Index t=findIndex(psi(RefBin-1),"bath"); //Print(t);
    ITensor Sm = ITensor(t,prime(t));
    Sm.set(t(2),prime(t(1)),1.);
    int counter=0;
    for(int i=RefBin;i>RefBin-SpectrumSteps;i--)
    {   
    t=findIndex(psi(i-1),"bath"); //Print(t);
    Sm = ITensor(t,prime(t));
    Sm.set(t(2),prime(t(1)),1.);

    temp = psi(i)*psi(i-1);

    Tg1= dag( noPrime(temp*Sp) )*noPrime(Sm*temp);
    Tg2= dag( noPrime(Sm*noPrime(Sp*temp) )) * noPrime(Sm*noPrime(Sp*temp));
    counter++;
    bath_g1[counter]=eltC(Tg1);
    bath_g2[counter]=eltC(Tg2);
    printf("bath_g1[%i]-|bc_inf|^2=%.10f+I%.10f -- g2=%.10f \n",counter,bath_g1[counter].real(),bath_g1[counter].imag(),bath_g2[counter].real());
    fprintf(file,"%.10f \t %.10f \t %10f \t %.10f \n",dt*counter,bath_g1[counter].real(),bath_g1[counter].imag(),bath_g2[counter].real());

    U = ITensor(bin[RefBin],commonIndex(psi(i-1),psi(i-2))); //PrintData(U);
    svd(temp,U,S,Tg1);
    psi.setA(i-1,U*S);
    psi.setA(i,Tg1); 
    }

     Complex b_offset = bath_g1[counter];
     double g2_norm = bath_g2[counter].real(); 
     

    fclose(file);    

// ---------------- output spectrum --------------------------------------
    
snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_TLS_Spectrum_w_FEEDBACK_at_Step=%i_of%i_OmL_%.3f_Gamma%.3f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,Nfb,(int)ttotal,Omega,Gamma,dt,cutoff*1000000.);
    file = fopen(FILE_NAME,"w");

    auto spectrum = std::vector<Complex>(SpectrumSteps+1);
   double dw = SpectrumIntervall/SpectrumSteps; 
   double om;
   
   for(int i=0;i<=SpectrumSteps;i++)
   {
     om = -0.5*(SpectrumIntervall)+i*dw;
     spectrum[i] = 0.;
     for(int b=0;b<=SpectrumSteps;b++) spectrum[i] += exp(-Cplx_i*om*b*dt)*(bath_g1[b]-b_offset);  
     fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \n",om,spectrum[i].real(),spectrum[i].imag(),bath_g2[i].real()/g2_norm);
   }

   fclose(file);


return 0;
}
