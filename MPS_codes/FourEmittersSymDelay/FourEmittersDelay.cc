#include "itensor/all.h"

using namespace itensor;


// ###########################################################################################################
// ############################### NO FB MPO #################################################################
// ###########################################################################################################
void MPO_NO_FB_SETUP
              (
               ITensor& U_evo, const std::vector<Index>& bin,
               int tls, int four_l, int Nbin,
               Real Gamma_l, Real Gamma_ml, Real Gamma_mr, Real Gamma_r, Real dt,
               Real phi_l, Real phi_ml,Real phi_mr, Real phi_r
              )
              
{
    
    int one_r = four_l+1;
    // left atom interaction hamiltonian ... time local field (now)
    // define: |abc> == Atom left in a, Atom middle b, Atom right c, Atom 3, Atom 2, Atom 1, 
    Complex   gam_l  =  -Cplx_i*Gamma_l*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_l);
    Complex c_gam_l  =  -Cplx_i*Gamma_l*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_l);

    Complex   gam_ml =  -Cplx_i*Gamma_ml*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_ml); 
    Complex c_gam_ml =  -Cplx_i*Gamma_ml*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_ml); 

    Complex   gam_mr =  -Cplx_i*Gamma_mr*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_mr); 
    Complex c_gam_mr =  -Cplx_i*Gamma_mr*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_mr); 

    Complex   gam_r  =  -Cplx_i*Gamma_r*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_r);
    Complex c_gam_r  =  -Cplx_i*Gamma_r*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_r);
   
    printf("INITIALIZE HAMILTONIAN: GAM_L=%.3f - GAM_ML=%.3f - GAM_MR=%.3f - GAM_R=%.3f\n",gam_l,gam_ml,gam_mr,gam_r);
    
    auto H_dis_l  = ITensor(bin[tls],prime(bin[tls]),bin[four_l],prime(bin[four_l]));
    auto H_dis_r  = ITensor(bin[tls],prime(bin[tls]),bin[one_r],prime(bin[one_r]));

    double sqrt_phot=1.;
    for( int phot=1;phot<Nbin;phot++)
    {    
    sqrt_phot=sqrt(1.*phot);    
    // -------------------------------- left atom -----------------------------------------------------------------
    // for the left atom deexcitation means -4 we go from unprimed to primed  
    //  |1><9|  + |2><10| + |3><11| + |5><13| + |7><15| + |6><14| + |4><12| + |8><16|       
    H_dis_l.set(bin[tls](9 ),prime(bin[tls](1 )),bin[four_l](phot),prime(bin[four_l](phot+1)),gam_l*sqrt_phot);
    H_dis_l.set(bin[tls](10),prime(bin[tls](2 )),bin[four_l](phot),prime(bin[four_l](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](11),prime(bin[tls](3 )),bin[four_l](phot),prime(bin[four_l](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](12),prime(bin[tls](4 )),bin[four_l](phot),prime(bin[four_l](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](13),prime(bin[tls](5 )),bin[four_l](phot),prime(bin[four_l](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](14),prime(bin[tls](6 )),bin[four_l](phot),prime(bin[four_l](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](15),prime(bin[tls](7 )),bin[four_l](phot),prime(bin[four_l](phot+1)),gam_l*sqrt_phot);
    H_dis_l.set(bin[tls](16),prime(bin[tls](8 )),bin[four_l](phot),prime(bin[four_l](phot+1)),gam_l*sqrt_phot);        
    // --------------------------------------------------------------------------------------------------------
    H_dis_l.set(bin[tls](1 ),prime(bin[tls](9 )),bin[four_l](phot+1),prime(bin[four_l](phot)),c_gam_l*sqrt_phot);
    H_dis_l.set(bin[tls](2 ),prime(bin[tls](10)),bin[four_l](phot+1),prime(bin[four_l](phot)),c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](3 ),prime(bin[tls](11)),bin[four_l](phot+1),prime(bin[four_l](phot)),c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](4 ),prime(bin[tls](12)),bin[four_l](phot+1),prime(bin[four_l](phot)),c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](5 ),prime(bin[tls](13)),bin[four_l](phot+1),prime(bin[four_l](phot)),c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](6 ),prime(bin[tls](14)),bin[four_l](phot+1),prime(bin[four_l](phot)),c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](7 ),prime(bin[tls](15)),bin[four_l](phot+1),prime(bin[four_l](phot)),c_gam_l*sqrt_phot);
    H_dis_l.set(bin[tls](8 ),prime(bin[tls](16)),bin[four_l](phot+1),prime(bin[four_l](phot)),c_gam_l*sqrt_phot);        
    // -------------------------------- middle left atom ------------------------------------------------------
    // |1><5|   + |2><6|   + |3><7| + |4><8| + |9><13| + |10><14| + |11><15|  + |12><16|
    H_dis_l.set(bin[tls](5 ),prime(bin[tls](1 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_ml*sqrt_phot);
    H_dis_l.set(bin[tls](6 ),prime(bin[tls](2 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_ml*sqrt_phot);        
    H_dis_l.set(bin[tls](7 ),prime(bin[tls](3 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_ml*sqrt_phot);        
    H_dis_l.set(bin[tls](8 ),prime(bin[tls](4 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_ml*sqrt_phot);        
    H_dis_l.set(bin[tls](13),prime(bin[tls](9 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_ml*sqrt_phot);
    H_dis_l.set(bin[tls](14),prime(bin[tls](10)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_ml*sqrt_phot);        
    H_dis_l.set(bin[tls](15),prime(bin[tls](11)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_ml*sqrt_phot);        
    H_dis_l.set(bin[tls](16),prime(bin[tls](12)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_ml*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_dis_l.set(bin[tls](1 ),prime(bin[tls](5 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_ml*sqrt_phot);
    H_dis_l.set(bin[tls](2 ),prime(bin[tls](6 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_ml*sqrt_phot);        
    H_dis_l.set(bin[tls](3 ),prime(bin[tls](7 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_ml*sqrt_phot);        
    H_dis_l.set(bin[tls](4 ),prime(bin[tls](8 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_ml*sqrt_phot);        
    H_dis_l.set(bin[tls](9 ),prime(bin[tls](13)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_ml*sqrt_phot);
    H_dis_l.set(bin[tls](10),prime(bin[tls](14)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_ml*sqrt_phot);        
    H_dis_l.set(bin[tls](11),prime(bin[tls](15)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_ml*sqrt_phot);        
    H_dis_l.set(bin[tls](12),prime(bin[tls](16)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_ml*sqrt_phot);            
    // -------------------------------- middle right atom -------------------------------------------------------
    // |1><3| + |2><4| + |5><7| + |6><8| + |9><11| + |10><12| + |13><15|  + |14><16|  
    H_dis_l.set(bin[tls](3 ),prime(bin[tls](1 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_mr*sqrt_phot);
    H_dis_l.set(bin[tls](4 ),prime(bin[tls](2 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_mr*sqrt_phot);        
    H_dis_l.set(bin[tls](7 ),prime(bin[tls](5 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_mr*sqrt_phot);        
    H_dis_l.set(bin[tls](8 ),prime(bin[tls](6 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_mr*sqrt_phot);        
    H_dis_l.set(bin[tls](11),prime(bin[tls](9 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_mr*sqrt_phot);
    H_dis_l.set(bin[tls](12),prime(bin[tls](10)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_mr*sqrt_phot);        
    H_dis_l.set(bin[tls](15),prime(bin[tls](13)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_mr*sqrt_phot);        
    H_dis_l.set(bin[tls](16),prime(bin[tls](14)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_mr*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_dis_l.set(bin[tls](1 ),prime(bin[tls](3 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_mr*sqrt_phot);
    H_dis_l.set(bin[tls](2 ),prime(bin[tls](4 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_mr*sqrt_phot);        
    H_dis_l.set(bin[tls](5 ),prime(bin[tls](7 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_mr*sqrt_phot);        
    H_dis_l.set(bin[tls](6 ),prime(bin[tls](8 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_mr*sqrt_phot);        
    H_dis_l.set(bin[tls](9 ),prime(bin[tls](11)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_mr*sqrt_phot);
    H_dis_l.set(bin[tls](10),prime(bin[tls](12)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_mr*sqrt_phot);        
    H_dis_l.set(bin[tls](13),prime(bin[tls](15)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_mr*sqrt_phot);        
    H_dis_l.set(bin[tls](14),prime(bin[tls](16)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_mr*sqrt_phot);            
    // -------------------------------- right atom -------------------------------------------------------
    // |1><2|   + |3><4|   + |5><6| + |7><8| + |9><10|+ |11><12| + |13><14| + |15><16|       
    H_dis_l.set(bin[tls](2 ),prime(bin[tls](1 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);
    H_dis_l.set(bin[tls](4 ),prime(bin[tls](3 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    H_dis_l.set(bin[tls](6 ),prime(bin[tls](5 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    H_dis_l.set(bin[tls](8 ),prime(bin[tls](7 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    H_dis_l.set(bin[tls](10),prime(bin[tls](9 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);
    H_dis_l.set(bin[tls](12),prime(bin[tls](11)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    H_dis_l.set(bin[tls](14),prime(bin[tls](13)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    H_dis_l.set(bin[tls](16),prime(bin[tls](15)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_dis_l.set(bin[tls](1 ),prime(bin[tls](2 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);
    H_dis_l.set(bin[tls](3 ),prime(bin[tls](4 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);        
    H_dis_l.set(bin[tls](5 ),prime(bin[tls](6 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);        
    H_dis_l.set(bin[tls](7 ),prime(bin[tls](8 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);        
    H_dis_l.set(bin[tls](9 ),prime(bin[tls](10)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);
    H_dis_l.set(bin[tls](11),prime(bin[tls](12)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);        
    H_dis_l.set(bin[tls](13),prime(bin[tls](14)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);        
    H_dis_l.set(bin[tls](15),prime(bin[tls](16)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);            

    // ################################ right moving field ###################################################### 

    // -------------------------------- left atom ---------------------------------------------------------------
    // for the left atom deexcitation means -4 we go from unprimed to primed  
    //  |1><9|  + |2><10| + |3><11| + |5><13| + |7><15| + |6><14| + |4><12| + |8><16|       
    H_dis_r.set(bin[tls](9 ),prime(bin[tls](1 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);
    H_dis_r.set(bin[tls](10),prime(bin[tls](2 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    H_dis_r.set(bin[tls](11),prime(bin[tls](3 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    H_dis_r.set(bin[tls](12),prime(bin[tls](4 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    H_dis_r.set(bin[tls](13),prime(bin[tls](5 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    H_dis_r.set(bin[tls](14),prime(bin[tls](6 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    H_dis_r.set(bin[tls](15),prime(bin[tls](7 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);
    H_dis_r.set(bin[tls](16),prime(bin[tls](8 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    // --------------------------------------------------------------------------------------------------------
    H_dis_r.set(bin[tls](1 ),prime(bin[tls](9 )),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);
    H_dis_r.set(bin[tls](2 ),prime(bin[tls](10)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    H_dis_r.set(bin[tls](3 ),prime(bin[tls](11)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    H_dis_r.set(bin[tls](4 ),prime(bin[tls](12)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    H_dis_r.set(bin[tls](5 ),prime(bin[tls](13)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    H_dis_r.set(bin[tls](6 ),prime(bin[tls](14)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    H_dis_r.set(bin[tls](7 ),prime(bin[tls](15)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);
    H_dis_r.set(bin[tls](8 ),prime(bin[tls](16)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    // -------------------------------- middle left atom ------------------------------------------------------
    // |1><5|   + |2><6|   + |3><7| + |4><8| + |9><13| + |10><14| + |11><15|  + |12><16|
    H_dis_r.set(bin[tls](5 ),prime(bin[tls](1 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_ml*sqrt_phot);
    H_dis_r.set(bin[tls](6 ),prime(bin[tls](2 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_ml*sqrt_phot);        
    H_dis_r.set(bin[tls](7 ),prime(bin[tls](3 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_ml*sqrt_phot);        
    H_dis_r.set(bin[tls](8 ),prime(bin[tls](4 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_ml*sqrt_phot);        
    H_dis_r.set(bin[tls](13),prime(bin[tls](9 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_ml*sqrt_phot);
    H_dis_r.set(bin[tls](14),prime(bin[tls](10)),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_ml*sqrt_phot);        
    H_dis_r.set(bin[tls](15),prime(bin[tls](11)),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_ml*sqrt_phot);        
    H_dis_r.set(bin[tls](16),prime(bin[tls](12)),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_ml*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_dis_r.set(bin[tls](1 ),prime(bin[tls](5 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_ml*sqrt_phot);
    H_dis_r.set(bin[tls](2 ),prime(bin[tls](6 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_ml*sqrt_phot);        
    H_dis_r.set(bin[tls](3 ),prime(bin[tls](7 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_ml*sqrt_phot);        
    H_dis_r.set(bin[tls](4 ),prime(bin[tls](8 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_ml*sqrt_phot);        
    H_dis_r.set(bin[tls](9 ),prime(bin[tls](13)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_ml*sqrt_phot);
    H_dis_r.set(bin[tls](10),prime(bin[tls](14)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_ml*sqrt_phot);        
    H_dis_r.set(bin[tls](11),prime(bin[tls](15)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_ml*sqrt_phot);        
    H_dis_r.set(bin[tls](12),prime(bin[tls](16)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_ml*sqrt_phot);            
    // -------------------------------- middle right atom -------------------------------------------------------
    // |1><3| + |2><4| + |5><7| + |6><8| + |9><11| + |10><12| + |13><15|  + |14><16|  
    H_dis_r.set(bin[tls](3 ),prime(bin[tls](1 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_mr*sqrt_phot);
    H_dis_r.set(bin[tls](4 ),prime(bin[tls](2 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_mr*sqrt_phot);        
    H_dis_r.set(bin[tls](7 ),prime(bin[tls](5 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_mr*sqrt_phot);        
    H_dis_r.set(bin[tls](8 ),prime(bin[tls](6 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_mr*sqrt_phot);        
    H_dis_r.set(bin[tls](11),prime(bin[tls](9 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_mr*sqrt_phot);
    H_dis_r.set(bin[tls](12),prime(bin[tls](10)),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_mr*sqrt_phot);        
    H_dis_r.set(bin[tls](15),prime(bin[tls](13)),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_mr*sqrt_phot);        
    H_dis_r.set(bin[tls](16),prime(bin[tls](14)),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_mr*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_dis_r.set(bin[tls](1 ),prime(bin[tls](3 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_mr*sqrt_phot);
    H_dis_r.set(bin[tls](2 ),prime(bin[tls](4 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_mr*sqrt_phot);        
    H_dis_r.set(bin[tls](5 ),prime(bin[tls](7 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_mr*sqrt_phot);        
    H_dis_r.set(bin[tls](6 ),prime(bin[tls](8 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_mr*sqrt_phot);        
    H_dis_r.set(bin[tls](9 ),prime(bin[tls](11)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_mr*sqrt_phot);
    H_dis_r.set(bin[tls](10),prime(bin[tls](12)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_mr*sqrt_phot);        
    H_dis_r.set(bin[tls](13),prime(bin[tls](15)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_mr*sqrt_phot);        
    H_dis_r.set(bin[tls](14),prime(bin[tls](16)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_mr*sqrt_phot);            
    // -------------------------------- right atom -------------------------------------------------------
    // |1><2|   + |3><4|   + |5><6| + |7><8| + |9><10|+ |11><12| + |13><14| + |15><16|       
    H_dis_r.set(bin[tls](2 ),prime(bin[tls](1 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_r*sqrt_phot);
    H_dis_r.set(bin[tls](4 ),prime(bin[tls](3 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](6 ),prime(bin[tls](5 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](8 ),prime(bin[tls](7 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](10),prime(bin[tls](9 )),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_r*sqrt_phot);
    H_dis_r.set(bin[tls](12),prime(bin[tls](11)),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](14),prime(bin[tls](13)),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](16),prime(bin[tls](15)),bin[one_r](phot),prime(bin[one_r](phot+1)), gam_r*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_dis_r.set(bin[tls](1 ),prime(bin[tls](2 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_r*sqrt_phot);
    H_dis_r.set(bin[tls](3 ),prime(bin[tls](4 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](5 ),prime(bin[tls](6 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](7 ),prime(bin[tls](8 )),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](9 ),prime(bin[tls](10)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_r*sqrt_phot);
    H_dis_r.set(bin[tls](11),prime(bin[tls](12)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](13),prime(bin[tls](14)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](15),prime(bin[tls](16)),bin[one_r](phot+1),prime(bin[one_r](phot)), c_gam_r*sqrt_phot);            
    }
    
    H_dis_r = delta(bin[four_l],prime(bin[four_l])) * H_dis_r;
    H_dis_l = delta(bin[one_r],prime(bin[one_r])) * H_dis_l;
    
    
    auto H_int = H_dis_r + H_dis_l ;  //PrintData(H_int);           
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

    Print(H_int);
    
    auto temp= ITensor(bin[tls],prime(bin[tls]));
    temp.set(bin[tls](1 ),prime(bin[tls](1 )),1.);
    temp.set(bin[tls](2 ),prime(bin[tls](2 )),1.);
    temp.set(bin[tls](3 ),prime(bin[tls](3 )),1.);
    temp.set(bin[tls](4 ),prime(bin[tls](4 )),1.);
    temp.set(bin[tls](5 ),prime(bin[tls](5 )),1.);
    temp.set(bin[tls](6 ),prime(bin[tls](6 )),1.);
    temp.set(bin[tls](7 ),prime(bin[tls](7 )),1.);
    temp.set(bin[tls](8 ),prime(bin[tls](8 )),1.);
    temp.set(bin[tls](9 ),prime(bin[tls](9 )),1.);
    temp.set(bin[tls](10),prime(bin[tls](10)),1.);
    temp.set(bin[tls](11),prime(bin[tls](11)),1.);
    temp.set(bin[tls](12),prime(bin[tls](12)),1.);
    temp.set(bin[tls](13),prime(bin[tls](13)),1.);
    temp.set(bin[tls](14),prime(bin[tls](14)),1.);
    temp.set(bin[tls](15),prime(bin[tls](15)),1.);
    temp.set(bin[tls](16),prime(bin[tls](16)),1.);

    //ITensor temp;
    temp = delta(bin[four_l],prime(bin[four_l]))*temp*delta(bin[one_r],prime(bin[one_r])); //Print(temp);  
    // now the taylor expansion of the U_evo

    U_evo =   temp + H_int_1 + H_int_2 + H_int_3 
                   + H_int_4 + H_int_5 + H_int_6 + H_int_7  + H_int_8 + H_int_9 + H_int_10 ;
                   
            
                   
    return;               
}
// ###########################################################################################################
// ###########################################################################################################
// ###########################################################################################################

// ###########################################################################################################
// ############################### FB MPO ####################################################################
// ###########################################################################################################
void MPO_SETUP(ITensor& U_evo, const std::vector<Index>& bin,
               int tls, int one_l, int two_l,int three_l, int four_l, int one_r, int two_r,int three_r, int four_r, int Nbin,
               Real Gamma_l, Real Gamma_ml, Real Gamma_mr, Real Gamma_r,Real phi_l, Real phi_ml, Real phi_mr, Real phi_r, Real dt
              )
{
    
    Complex   gam_l    =  -Cplx_i*Gamma_l*sqrt(dt);
    Complex c_gam_l    =  -Cplx_i*Gamma_l*sqrt(dt);
    Complex   gam_l_fb =  -Cplx_i*Gamma_l*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_l);
    Complex c_gam_l_fb =  -Cplx_i*Gamma_l*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_l);

    Complex   gam_mll_fb =  -Cplx_i*Gamma_ml*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_ml); 
    Complex c_gam_mll_fb =  -Cplx_i*Gamma_ml*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_ml);  
    Complex   gam_mlr_fb =  -Cplx_i*Gamma_ml*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_mr); 
    Complex c_gam_mlr_fb =  -Cplx_i*Gamma_ml*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_mr);  

    // third atom, middle-left picks up phi*1./3. with left moving field
    Complex   gam_mrl_fb =  -Cplx_i*Gamma_mr*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_mr); 
    Complex c_gam_mrl_fb =  -Cplx_i*Gamma_mr*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_mr); 
    // third atom, middle-right, picks up phi*2./3. with right moving field
    Complex   gam_mrr_fb =  -Cplx_i*Gamma_mr*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_ml); 
    Complex c_gam_mrr_fb =  -Cplx_i*Gamma_mr*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_ml);   

    Complex   gam_r    =  -Cplx_i*Gamma_r*sqrt(dt);
    Complex c_gam_r    =  -Cplx_i*Gamma_r*sqrt(dt);  
    Complex   gam_r_fb =  -Cplx_i*Gamma_r*sqrt(dt)*exp(-Cplx_i*3.14159265359*phi_r); 
    Complex c_gam_r_fb =  -Cplx_i*Gamma_r*sqrt(dt)*exp( Cplx_i*3.14159265359*phi_r);
   
    printf("INITIALIZE HAMILTONIAN: GAM_L  =%.3f - GAM_L_FB=%.3f - GAM_R  =%.3f - GAM_R_FB=%.3f\n",gam_l,gam_l_fb,gam_r,gam_r_fb);
    printf("INITIALIZE HAMILTONIAN: GAM_MLL=%.3f - GAM_MLR =%.3f - GAM_MRL=%.3f - GAM_RR  =%.3f\n",gam_mll_fb,gam_mlr_fb,gam_mrl_fb,gam_mrr_fb);
 
    
    // two-time local interaction hamiltonians for emitter at left and right border
    auto H_dis_l  = ITensor(bin[tls],prime(bin[tls]),bin[one_r],prime(bin[one_r]));
    auto H_dis_r  = ITensor(bin[tls],prime(bin[tls]),bin[four_l],prime(bin[four_l]));
    
    double sqrt_phot=1.;
    for( int phot=1;phot<Nbin;phot++)
    {    
    sqrt_phot=sqrt(1.*phot);    
    // -------------------------------- left atom -----------------------------------------------------------------
    // for the left atom deexcitation means -4 we go from unprimed to primed  
    //  |1><9|  + |2><10| + |3><11| + |5><13| + |7><15| + |6><14| + |4><12| + |8><16|       
    H_dis_l.set(bin[tls](9 ),prime(bin[tls](1 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);
    H_dis_l.set(bin[tls](10),prime(bin[tls](2 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](11),prime(bin[tls](3 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](12),prime(bin[tls](4 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](13),prime(bin[tls](5 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](14),prime(bin[tls](6 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](15),prime(bin[tls](7 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);
    H_dis_l.set(bin[tls](16),prime(bin[tls](8 )),bin[one_r](phot),prime(bin[one_r](phot+1)),gam_l*sqrt_phot);        
    // --------------------------------------------------------------------------------------------------------
    H_dis_l.set(bin[tls](1 ),prime(bin[tls](9 )),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);
    H_dis_l.set(bin[tls](2 ),prime(bin[tls](10)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](3 ),prime(bin[tls](11)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](4 ),prime(bin[tls](12)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](5 ),prime(bin[tls](13)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](6 ),prime(bin[tls](14)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    H_dis_l.set(bin[tls](7 ),prime(bin[tls](15)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);
    H_dis_l.set(bin[tls](8 ),prime(bin[tls](16)),bin[one_r](phot+1),prime(bin[one_r](phot)),c_gam_l*sqrt_phot);        
    // -------------------------------- right atom ------------------------------------------------------------
    // ----------- |1><2|   + |3><4|   + |5><6| + |7><8| + |9><10|+ |11><12| + |13><14| + |15><16|   ----------    
    // -------------------------------- right atom ------------------------------------------------------------
    H_dis_r.set(bin[tls](2 ),prime(bin[tls](1 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);
    H_dis_r.set(bin[tls](4 ),prime(bin[tls](3 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](6 ),prime(bin[tls](5 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](8 ),prime(bin[tls](7 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](10),prime(bin[tls](9 )),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);
    H_dis_r.set(bin[tls](12),prime(bin[tls](11)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](14),prime(bin[tls](13)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](16),prime(bin[tls](15)),bin[four_l](phot),prime(bin[four_l](phot+1)), gam_r*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_dis_r.set(bin[tls](1 ),prime(bin[tls](2 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);
    H_dis_r.set(bin[tls](3 ),prime(bin[tls](4 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](5 ),prime(bin[tls](6 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](7 ),prime(bin[tls](8 )),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](9 ),prime(bin[tls](10)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);
    H_dis_r.set(bin[tls](11),prime(bin[tls](12)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](13),prime(bin[tls](14)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot);        
    H_dis_r.set(bin[tls](15),prime(bin[tls](16)),bin[four_l](phot+1),prime(bin[four_l](phot)), c_gam_r*sqrt_phot); 
    }
    
    // now bring dissipation hamiltonians into one product Hilbert space
    H_dis_l = delta(bin[four_l     ],prime(bin[four_l     ])) * H_dis_l;
    H_dis_l = delta(bin[two_r],prime(bin[two_r])) * H_dis_l;
    H_dis_l = delta(bin[two_l],prime(bin[two_l])) * H_dis_l;
    H_dis_l = delta(bin[three_r],prime(bin[three_r])) * H_dis_l;
    H_dis_l = delta(bin[three_l],prime(bin[three_l])) * H_dis_l;
    H_dis_l = delta(bin[one_l    ],prime(bin[one_l    ])) * H_dis_l;
    H_dis_l = delta(bin[four_r    ],prime(bin[four_r    ])) * H_dis_l;

    H_dis_r = delta(bin[one_r     ],prime(bin[one_r     ])) * H_dis_r;
    H_dis_r = delta(bin[two_r],prime(bin[two_r])) * H_dis_r;
    H_dis_r = delta(bin[two_l],prime(bin[two_l])) * H_dis_r;
    H_dis_r = delta(bin[three_r],prime(bin[three_r])) * H_dis_r;
    H_dis_r = delta(bin[three_l],prime(bin[three_l])) * H_dis_r;
    H_dis_r = delta(bin[one_l    ],prime(bin[one_l    ])) * H_dis_r;
    H_dis_r = delta(bin[four_r    ],prime(bin[four_r    ])) * H_dis_r;
    
    // Print(H_dis_l);Print(H_dis_r);
    
    // time non-local field with excitation exchange for left and right atom as the middle one has no extra feedback 
    auto H_ml_l = ITensor(bin[tls],prime(bin[tls]),bin[two_l],prime(bin[two_l]));
    auto H_ml_r = ITensor(bin[tls],prime(bin[tls]),bin[two_r],prime(bin[two_r]));

    auto H_mr_l = ITensor(bin[tls],prime(bin[tls]),bin[three_l],prime(bin[three_l]));
    auto H_mr_r = ITensor(bin[tls],prime(bin[tls]),bin[three_r],prime(bin[three_r]));

    auto H_fb_l  = ITensor(bin[tls],prime(bin[tls]),bin[one_l],prime(bin[one_l]));
    auto H_fb_r  = ITensor(bin[tls],prime(bin[tls]),bin[four_r],prime(bin[four_r]));
    
    for( int phot=1;phot<Nbin;phot++)
    {   
    sqrt_phot=sqrt(1.*phot);        
    // for the left atom deexcitation means -4 we go from unprimed to primed  
    //  |1><9|  + |2><10| + |3><11| + |5><13| + |7><15| + |6><14| + |4><12| + |8><16|       
    H_fb_l.set(bin[tls](9 ),prime(bin[tls](1 )),bin[one_l](phot),prime(bin[one_l](phot+1)),gam_l_fb*sqrt_phot);
    H_fb_l.set(bin[tls](10),prime(bin[tls](2 )),bin[one_l](phot),prime(bin[one_l](phot+1)),gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](11),prime(bin[tls](3 )),bin[one_l](phot),prime(bin[one_l](phot+1)),gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](12),prime(bin[tls](4 )),bin[one_l](phot),prime(bin[one_l](phot+1)),gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](13),prime(bin[tls](5 )),bin[one_l](phot),prime(bin[one_l](phot+1)),gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](14),prime(bin[tls](6 )),bin[one_l](phot),prime(bin[one_l](phot+1)),gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](15),prime(bin[tls](7 )),bin[one_l](phot),prime(bin[one_l](phot+1)),gam_l_fb*sqrt_phot);
    H_fb_l.set(bin[tls](16),prime(bin[tls](8 )),bin[one_l](phot),prime(bin[one_l](phot+1)),gam_l_fb*sqrt_phot);        
    // --------------------------------------------------------------------------------------------------------
    H_fb_l.set(bin[tls](1 ),prime(bin[tls](9 )),bin[one_l](phot+1),prime(bin[one_l](phot)),c_gam_l_fb*sqrt_phot);
    H_fb_l.set(bin[tls](2 ),prime(bin[tls](10)),bin[one_l](phot+1),prime(bin[one_l](phot)),c_gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](3 ),prime(bin[tls](11)),bin[one_l](phot+1),prime(bin[one_l](phot)),c_gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](4 ),prime(bin[tls](12)),bin[one_l](phot+1),prime(bin[one_l](phot)),c_gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](5 ),prime(bin[tls](13)),bin[one_l](phot+1),prime(bin[one_l](phot)),c_gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](6 ),prime(bin[tls](14)),bin[one_l](phot+1),prime(bin[one_l](phot)),c_gam_l_fb*sqrt_phot);        
    H_fb_l.set(bin[tls](7 ),prime(bin[tls](15)),bin[one_l](phot+1),prime(bin[one_l](phot)),c_gam_l_fb*sqrt_phot);
    H_fb_l.set(bin[tls](8 ),prime(bin[tls](16)),bin[one_l](phot+1),prime(bin[one_l](phot)),c_gam_l_fb*sqrt_phot);  
    // -------------------------------- middle left atom ------------------------------------------------------
    // ---- |1><5|   + |2><6|   + |3><7| + |4><8| + |9><13| + |10><14| + |11><15|  + |12><16| -----------------
    // ------------------- left moving field ------------------------------------------------------------------
    H_ml_l.set(bin[tls](5 ),prime(bin[tls](1 )),bin[two_l](phot),prime(bin[two_l](phot+1)), gam_mll_fb*sqrt_phot);
    H_ml_l.set(bin[tls](6 ),prime(bin[tls](2 )),bin[two_l](phot),prime(bin[two_l](phot+1)), gam_mll_fb*sqrt_phot);        
    H_ml_l.set(bin[tls](7 ),prime(bin[tls](3 )),bin[two_l](phot),prime(bin[two_l](phot+1)), gam_mll_fb*sqrt_phot);        
    H_ml_l.set(bin[tls](8 ),prime(bin[tls](4 )),bin[two_l](phot),prime(bin[two_l](phot+1)), gam_mll_fb*sqrt_phot);        
    H_ml_l.set(bin[tls](13),prime(bin[tls](9 )),bin[two_l](phot),prime(bin[two_l](phot+1)), gam_mll_fb*sqrt_phot);
    H_ml_l.set(bin[tls](14),prime(bin[tls](10)),bin[two_l](phot),prime(bin[two_l](phot+1)), gam_mll_fb*sqrt_phot);        
    H_ml_l.set(bin[tls](15),prime(bin[tls](11)),bin[two_l](phot),prime(bin[two_l](phot+1)), gam_mll_fb*sqrt_phot);        
    H_ml_l.set(bin[tls](16),prime(bin[tls](12)),bin[two_l](phot),prime(bin[two_l](phot+1)), gam_mll_fb*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_ml_l.set(bin[tls](1 ),prime(bin[tls](5 )),bin[two_l](phot+1),prime(bin[two_l](phot)), c_gam_mll_fb*sqrt_phot);
    H_ml_l.set(bin[tls](2 ),prime(bin[tls](6 )),bin[two_l](phot+1),prime(bin[two_l](phot)), c_gam_mll_fb*sqrt_phot);        
    H_ml_l.set(bin[tls](3 ),prime(bin[tls](7 )),bin[two_l](phot+1),prime(bin[two_l](phot)), c_gam_mll_fb*sqrt_phot);        
    H_ml_l.set(bin[tls](4 ),prime(bin[tls](8 )),bin[two_l](phot+1),prime(bin[two_l](phot)), c_gam_mll_fb*sqrt_phot);        
    H_ml_l.set(bin[tls](9 ),prime(bin[tls](13)),bin[two_l](phot+1),prime(bin[two_l](phot)), c_gam_mll_fb*sqrt_phot);
    H_ml_l.set(bin[tls](10),prime(bin[tls](14)),bin[two_l](phot+1),prime(bin[two_l](phot)), c_gam_mll_fb*sqrt_phot);        
    H_ml_l.set(bin[tls](11),prime(bin[tls](15)),bin[two_l](phot+1),prime(bin[two_l](phot)), c_gam_mll_fb*sqrt_phot);        
    H_ml_l.set(bin[tls](12),prime(bin[tls](16)),bin[two_l](phot+1),prime(bin[two_l](phot)), c_gam_mll_fb*sqrt_phot);            
    // ------------------- right moving field ------------------------------------------------------------------
    H_ml_r.set(bin[tls](5 ),prime(bin[tls](1 )),bin[two_r](phot),prime(bin[two_r](phot+1)), gam_mlr_fb*sqrt_phot);
    H_ml_r.set(bin[tls](6 ),prime(bin[tls](2 )),bin[two_r](phot),prime(bin[two_r](phot+1)), gam_mlr_fb*sqrt_phot);        
    H_ml_r.set(bin[tls](7 ),prime(bin[tls](3 )),bin[two_r](phot),prime(bin[two_r](phot+1)), gam_mlr_fb*sqrt_phot);        
    H_ml_r.set(bin[tls](8 ),prime(bin[tls](4 )),bin[two_r](phot),prime(bin[two_r](phot+1)), gam_mlr_fb*sqrt_phot);        
    H_ml_r.set(bin[tls](13),prime(bin[tls](9 )),bin[two_r](phot),prime(bin[two_r](phot+1)), gam_mlr_fb*sqrt_phot);
    H_ml_r.set(bin[tls](14),prime(bin[tls](10)),bin[two_r](phot),prime(bin[two_r](phot+1)), gam_mlr_fb*sqrt_phot);        
    H_ml_r.set(bin[tls](15),prime(bin[tls](11)),bin[two_r](phot),prime(bin[two_r](phot+1)), gam_mlr_fb*sqrt_phot);        
    H_ml_r.set(bin[tls](16),prime(bin[tls](12)),bin[two_r](phot),prime(bin[two_r](phot+1)), gam_mlr_fb*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_ml_r.set(bin[tls](1 ),prime(bin[tls](5 )),bin[two_r](phot+1),prime(bin[two_r](phot)), c_gam_mlr_fb*sqrt_phot);
    H_ml_r.set(bin[tls](2 ),prime(bin[tls](6 )),bin[two_r](phot+1),prime(bin[two_r](phot)), c_gam_mlr_fb*sqrt_phot);        
    H_ml_r.set(bin[tls](3 ),prime(bin[tls](7 )),bin[two_r](phot+1),prime(bin[two_r](phot)), c_gam_mlr_fb*sqrt_phot);        
    H_ml_r.set(bin[tls](4 ),prime(bin[tls](8 )),bin[two_r](phot+1),prime(bin[two_r](phot)), c_gam_mlr_fb*sqrt_phot);        
    H_ml_r.set(bin[tls](9 ),prime(bin[tls](13)),bin[two_r](phot+1),prime(bin[two_r](phot)), c_gam_mlr_fb*sqrt_phot);
    H_ml_r.set(bin[tls](10),prime(bin[tls](14)),bin[two_r](phot+1),prime(bin[two_r](phot)), c_gam_mlr_fb*sqrt_phot);        
    H_ml_r.set(bin[tls](11),prime(bin[tls](15)),bin[two_r](phot+1),prime(bin[two_r](phot)), c_gam_mlr_fb*sqrt_phot);        
    H_ml_r.set(bin[tls](12),prime(bin[tls](16)),bin[two_r](phot+1),prime(bin[two_r](phot)), c_gam_mlr_fb*sqrt_phot);            
    // -------------------------------- middle right atom ------------------------------------------------------
    // ---- |1><5|   + |2><6|   + |3><7| + |4><8| + |9><13| + |10><14| + |11><15|  + |12><16| -----------------
    // ------------------- left moving field ------------------------------------------------------------------
    H_mr_l.set(bin[tls](3 ),prime(bin[tls](1 )),bin[three_l](phot),prime(bin[three_l](phot+1)), gam_mrl_fb*sqrt_phot);
    H_mr_l.set(bin[tls](4 ),prime(bin[tls](2 )),bin[three_l](phot),prime(bin[three_l](phot+1)), gam_mrl_fb*sqrt_phot);        
    H_mr_l.set(bin[tls](7 ),prime(bin[tls](5 )),bin[three_l](phot),prime(bin[three_l](phot+1)), gam_mrl_fb*sqrt_phot);        
    H_mr_l.set(bin[tls](8 ),prime(bin[tls](6 )),bin[three_l](phot),prime(bin[three_l](phot+1)), gam_mrl_fb*sqrt_phot);        
    H_mr_l.set(bin[tls](11),prime(bin[tls](9 )),bin[three_l](phot),prime(bin[three_l](phot+1)), gam_mrl_fb*sqrt_phot);
    H_mr_l.set(bin[tls](12),prime(bin[tls](10)),bin[three_l](phot),prime(bin[three_l](phot+1)), gam_mrl_fb*sqrt_phot);        
    H_mr_l.set(bin[tls](15),prime(bin[tls](13)),bin[three_l](phot),prime(bin[three_l](phot+1)), gam_mrl_fb*sqrt_phot);        
    H_mr_l.set(bin[tls](16),prime(bin[tls](14)),bin[three_l](phot),prime(bin[three_l](phot+1)), gam_mrl_fb*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_mr_l.set(bin[tls](1 ),prime(bin[tls](3 )),bin[three_l](phot+1),prime(bin[three_l](phot)), c_gam_mrl_fb*sqrt_phot);
    H_mr_l.set(bin[tls](2 ),prime(bin[tls](4 )),bin[three_l](phot+1),prime(bin[three_l](phot)), c_gam_mrl_fb*sqrt_phot);        
    H_mr_l.set(bin[tls](5 ),prime(bin[tls](7 )),bin[three_l](phot+1),prime(bin[three_l](phot)), c_gam_mrl_fb*sqrt_phot);        
    H_mr_l.set(bin[tls](6 ),prime(bin[tls](8 )),bin[three_l](phot+1),prime(bin[three_l](phot)), c_gam_mrl_fb*sqrt_phot);        
    H_mr_l.set(bin[tls](9 ),prime(bin[tls](11)),bin[three_l](phot+1),prime(bin[three_l](phot)), c_gam_mrl_fb*sqrt_phot);
    H_mr_l.set(bin[tls](10),prime(bin[tls](12)),bin[three_l](phot+1),prime(bin[three_l](phot)), c_gam_mrl_fb*sqrt_phot);        
    H_mr_l.set(bin[tls](13),prime(bin[tls](15)),bin[three_l](phot+1),prime(bin[three_l](phot)), c_gam_mrl_fb*sqrt_phot);        
    H_mr_l.set(bin[tls](14),prime(bin[tls](16)),bin[three_l](phot+1),prime(bin[three_l](phot)), c_gam_mrl_fb*sqrt_phot);            
    // -------------------------------- right moving field ------------------------------------------------------------------
    H_mr_r.set(bin[tls](3 ),prime(bin[tls](1 )),bin[three_r](phot),prime(bin[three_r](phot+1)), gam_mrr_fb*sqrt_phot);
    H_mr_r.set(bin[tls](4 ),prime(bin[tls](2 )),bin[three_r](phot),prime(bin[three_r](phot+1)), gam_mrr_fb*sqrt_phot);        
    H_mr_r.set(bin[tls](7 ),prime(bin[tls](5 )),bin[three_r](phot),prime(bin[three_r](phot+1)), gam_mrr_fb*sqrt_phot);        
    H_mr_r.set(bin[tls](8 ),prime(bin[tls](6 )),bin[three_r](phot),prime(bin[three_r](phot+1)), gam_mrr_fb*sqrt_phot);        
    H_mr_r.set(bin[tls](11),prime(bin[tls](9 )),bin[three_r](phot),prime(bin[three_r](phot+1)), gam_mrr_fb*sqrt_phot);
    H_mr_r.set(bin[tls](12),prime(bin[tls](10)),bin[three_r](phot),prime(bin[three_r](phot+1)), gam_mrr_fb*sqrt_phot);        
    H_mr_r.set(bin[tls](15),prime(bin[tls](13)),bin[three_r](phot),prime(bin[three_r](phot+1)), gam_mrr_fb*sqrt_phot);        
    H_mr_r.set(bin[tls](16),prime(bin[tls](14)),bin[three_r](phot),prime(bin[three_r](phot+1)), gam_mrr_fb*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_mr_r.set(bin[tls](1 ),prime(bin[tls](3 )),bin[three_r](phot+1),prime(bin[three_r](phot)), c_gam_mrr_fb*sqrt_phot);
    H_mr_r.set(bin[tls](2 ),prime(bin[tls](4 )),bin[three_r](phot+1),prime(bin[three_r](phot)), c_gam_mrr_fb*sqrt_phot);        
    H_mr_r.set(bin[tls](5 ),prime(bin[tls](7 )),bin[three_r](phot+1),prime(bin[three_r](phot)), c_gam_mrr_fb*sqrt_phot);        
    H_mr_r.set(bin[tls](6 ),prime(bin[tls](8 )),bin[three_r](phot+1),prime(bin[three_r](phot)), c_gam_mrr_fb*sqrt_phot);        
    H_mr_r.set(bin[tls](9 ),prime(bin[tls](11)),bin[three_r](phot+1),prime(bin[three_r](phot)), c_gam_mrr_fb*sqrt_phot);
    H_mr_r.set(bin[tls](10),prime(bin[tls](12)),bin[three_r](phot+1),prime(bin[three_r](phot)), c_gam_mrr_fb*sqrt_phot);        
    H_mr_r.set(bin[tls](13),prime(bin[tls](15)),bin[three_r](phot+1),prime(bin[three_r](phot)), c_gam_mrr_fb*sqrt_phot);        
    H_mr_r.set(bin[tls](14),prime(bin[tls](16)),bin[three_r](phot+1),prime(bin[three_r](phot)), c_gam_mrr_fb*sqrt_phot);
    // -------------------------------- right atom ------------------------------------------------------------
    // ----------- |1><2|   + |3><4|   + |5><6| + |7><8| + |9><10|+ |11><12| + |13><14| + |15><16|   ----------    
    // -------------------------------- right atom ------------------------------------------------------------
    H_fb_r.set(bin[tls](2 ),prime(bin[tls](1 )),bin[four_r](phot),prime(bin[four_r](phot+1)), gam_r_fb*sqrt_phot);
    H_fb_r.set(bin[tls](4 ),prime(bin[tls](3 )),bin[four_r](phot),prime(bin[four_r](phot+1)), gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](6 ),prime(bin[tls](5 )),bin[four_r](phot),prime(bin[four_r](phot+1)), gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](8 ),prime(bin[tls](7 )),bin[four_r](phot),prime(bin[four_r](phot+1)), gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](10),prime(bin[tls](9 )),bin[four_r](phot),prime(bin[four_r](phot+1)), gam_r_fb*sqrt_phot);
    H_fb_r.set(bin[tls](12),prime(bin[tls](11)),bin[four_r](phot),prime(bin[four_r](phot+1)), gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](14),prime(bin[tls](13)),bin[four_r](phot),prime(bin[four_r](phot+1)), gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](16),prime(bin[tls](15)),bin[four_r](phot),prime(bin[four_r](phot+1)), gam_r_fb*sqrt_phot);        
    // ----------------------------------------------------------------------------------------------------------
    H_fb_r.set(bin[tls](1 ),prime(bin[tls](2 )),bin[four_r](phot+1),prime(bin[four_r](phot)), c_gam_r_fb*sqrt_phot);
    H_fb_r.set(bin[tls](3 ),prime(bin[tls](4 )),bin[four_r](phot+1),prime(bin[four_r](phot)), c_gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](5 ),prime(bin[tls](6 )),bin[four_r](phot+1),prime(bin[four_r](phot)), c_gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](7 ),prime(bin[tls](8 )),bin[four_r](phot+1),prime(bin[four_r](phot)), c_gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](9 ),prime(bin[tls](10)),bin[four_r](phot+1),prime(bin[four_r](phot)), c_gam_r_fb*sqrt_phot);
    H_fb_r.set(bin[tls](11),prime(bin[tls](12)),bin[four_r](phot+1),prime(bin[four_r](phot)), c_gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](13),prime(bin[tls](14)),bin[four_r](phot+1),prime(bin[four_r](phot)), c_gam_r_fb*sqrt_phot);        
    H_fb_r.set(bin[tls](15),prime(bin[tls](16)),bin[four_r](phot+1),prime(bin[four_r](phot)), c_gam_r_fb*sqrt_phot); 
    }
    
    H_ml_l = delta(bin[two_r],prime(bin[two_r])) * H_ml_l                                                ;
    H_ml_l = delta(bin[one_r     ],prime(bin[one_r     ])) * H_ml_l * delta(bin[four_l     ],prime(bin[four_l     ]));
    H_ml_l = delta(bin[three_r],prime(bin[three_r])) * H_ml_l * delta(bin[three_l],prime(bin[three_l]));
    H_ml_l = delta(bin[four_r    ],prime(bin[four_r    ])) * H_ml_l * delta(bin[one_l    ],prime(bin[one_l    ]));

    H_ml_r = delta(bin[two_l],prime(bin[two_l])) * H_ml_r                                                ;
    H_ml_r = delta(bin[one_r     ],prime(bin[one_r     ])) * H_ml_r * delta(bin[four_l     ],prime(bin[four_l     ]));
    H_ml_r = delta(bin[three_r],prime(bin[three_r])) * H_ml_r * delta(bin[three_l],prime(bin[three_l]));
    H_ml_r = delta(bin[four_r    ],prime(bin[four_r    ])) * H_ml_r * delta(bin[one_l    ],prime(bin[one_l    ]));

    H_mr_l = delta(bin[three_r],prime(bin[three_r])) * H_mr_l                                                ;
    H_mr_l = delta(bin[one_r     ],prime(bin[one_r     ])) * H_mr_l * delta(bin[four_l     ],prime(bin[four_l     ]));
    H_mr_l = delta(bin[two_r],prime(bin[two_r])) * H_mr_l * delta(bin[two_l],prime(bin[two_l]));
    H_mr_l = delta(bin[four_r    ],prime(bin[four_r    ])) * H_mr_l * delta(bin[one_l    ],prime(bin[one_l    ]));

    H_mr_r = delta(bin[three_l],prime(bin[three_l])) * H_mr_r                                                ;
    H_mr_r = delta(bin[one_r     ],prime(bin[one_r     ])) * H_mr_r * delta(bin[four_l     ],prime(bin[four_l     ]));
    H_mr_r = delta(bin[two_r],prime(bin[two_r])) * H_mr_r * delta(bin[two_l],prime(bin[two_l]));
    H_mr_r = delta(bin[four_r    ],prime(bin[four_r    ])) * H_mr_r * delta(bin[one_l    ],prime(bin[one_l    ]));

    H_fb_l = delta(bin[four_r    ],prime(bin[four_r    ])) * H_fb_l                                                ;
    H_fb_l = delta(bin[four_l     ],prime(bin[four_l     ])) * H_fb_l * delta(bin[one_r     ],prime(bin[one_r     ])); 
    H_fb_l = delta(bin[two_l],prime(bin[two_l])) * H_fb_l * delta(bin[two_r],prime(bin[two_r])); 
    H_fb_l = delta(bin[three_l],prime(bin[three_l])) * H_fb_l * delta(bin[three_r],prime(bin[three_r])); 

    H_fb_r = delta(bin[one_l    ],prime(bin[one_l    ])) * H_fb_r                                                ;
    H_fb_r = delta(bin[four_l     ],prime(bin[four_l     ])) * H_fb_r * delta(bin[one_r     ],prime(bin[one_r     ])); 
    H_fb_r = delta(bin[two_l],prime(bin[two_l])) * H_fb_r * delta(bin[two_r],prime(bin[two_r])); 
    H_fb_r = delta(bin[three_l],prime(bin[three_l])) * H_fb_r * delta(bin[three_r],prime(bin[three_r])); 
              
    //Print(H_dis_l);Print(H_dis_r);Print(H_fb_l);Print(H_fb_r); Print(H_ml_l);Print(H_ml_r); Print(H_mr_l);Print(H_mr_r);
    
    auto H_int = H_dis_l + H_dis_r + H_fb_l + H_fb_r + H_ml_l + H_ml_r + H_mr_l + H_mr_r ;  //PrintData(H_int);           
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
    auto H_int_10 = (1./10.) * mapPrime(temp_H_int_1*H_int_9,2,1); printf("10th-order prepared: order of H=%i.\n",order(H_int_10)); //*/

    auto temp= ITensor(bin[tls],prime(bin[tls]));
    temp.set(bin[tls](1 ),prime(bin[tls](1 )),1.);
    temp.set(bin[tls](2 ),prime(bin[tls](2 )),1.);
    temp.set(bin[tls](3 ),prime(bin[tls](3 )),1.);
    temp.set(bin[tls](4 ),prime(bin[tls](4 )),1.);
    temp.set(bin[tls](5 ),prime(bin[tls](5 )),1.);
    temp.set(bin[tls](6 ),prime(bin[tls](6 )),1.);
    temp.set(bin[tls](7 ),prime(bin[tls](7 )),1.);
    temp.set(bin[tls](8 ),prime(bin[tls](8 )),1.);
    temp.set(bin[tls](9 ),prime(bin[tls](9 )),1.);
    temp.set(bin[tls](10),prime(bin[tls](10)),1.);
    temp.set(bin[tls](11),prime(bin[tls](11)),1.);
    temp.set(bin[tls](12),prime(bin[tls](12)),1.);
    temp.set(bin[tls](13),prime(bin[tls](13)),1.);
    temp.set(bin[tls](14),prime(bin[tls](14)),1.);
    temp.set(bin[tls](15),prime(bin[tls](15)),1.);
    temp.set(bin[tls](16),prime(bin[tls](16)),1.);

    //ITensor temp;
    temp = delta(bin[four_l ],prime(bin[four_l ]))*temp*delta(bin[one_l  ],prime(bin[one_l  ]));
    temp = delta(bin[one_r  ],prime(bin[one_r  ]))*temp*delta(bin[four_r ],prime(bin[four_r ]));
    temp = delta(bin[two_l  ],prime(bin[two_l  ]))*temp*delta(bin[two_r  ],prime(bin[two_r  ]));
    temp = delta(bin[three_l],prime(bin[three_l]))*temp*delta(bin[three_r],prime(bin[three_r]));

    // now the taylor expansion of the U_evo
    U_evo =   temp + H_int_1 + H_int_2 + H_int_3 + H_int_4 + H_int_5 + H_int_6
                   + H_int_7  + H_int_8 + H_int_9 + H_int_10 
                   ;
                   
    return;               
}
// ###########################################################################################################
// ###########################################################################################################
// ###########################################################################################################
//*/
// ###########################################################################################################
// ###########################################################################################################
// ###########################################################################################################
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
    MPSTensor.set(bin[tls](1 ),binlink[Nfb](1),binlink[Nfb+1](1),Init[1 ]); // 0000 
    MPSTensor.set(bin[tls](2 ),binlink[Nfb](1),binlink[Nfb+1](1),Init[2 ]); // 0001
    MPSTensor.set(bin[tls](3 ),binlink[Nfb](1),binlink[Nfb+1](1),Init[3 ]); // 0010
    MPSTensor.set(bin[tls](4 ),binlink[Nfb](1),binlink[Nfb+1](1),Init[4 ]); // 0011
    MPSTensor.set(bin[tls](5 ),binlink[Nfb](1),binlink[Nfb+1](1),Init[5 ]); // 0100
    MPSTensor.set(bin[tls](6 ),binlink[Nfb](1),binlink[Nfb+1](1),Init[6 ]); // 0101
    MPSTensor.set(bin[tls](7 ),binlink[Nfb](1),binlink[Nfb+1](1),Init[7 ]); // 0110
    MPSTensor.set(bin[tls](8 ),binlink[Nfb](1),binlink[Nfb+1](1),Init[8 ]); // 0111
    MPSTensor.set(bin[tls](9 ),binlink[Nfb](1),binlink[Nfb+1](1),Init[9 ]); // 1000 
    MPSTensor.set(bin[tls](10),binlink[Nfb](1),binlink[Nfb+1](1),Init[10]); // 1001
    MPSTensor.set(bin[tls](11),binlink[Nfb](1),binlink[Nfb+1](1),Init[11]); // 1010
    MPSTensor.set(bin[tls](12),binlink[Nfb](1),binlink[Nfb+1](1),Init[12]); // 1011
    MPSTensor.set(bin[tls](13),binlink[Nfb](1),binlink[Nfb+1](1),Init[13]); // 1100
    MPSTensor.set(bin[tls](14),binlink[Nfb](1),binlink[Nfb+1](1),Init[14]); // 1101
    MPSTensor.set(bin[tls](15),binlink[Nfb](1),binlink[Nfb+1](1),Init[15]); // 1110
    MPSTensor.set(bin[tls](16),binlink[Nfb](1),binlink[Nfb+1](1),Init[16]); // 1111
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
// ###########################################################################################################
// ###########################################################################################################
// ###########################################################################################################

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

// ###########################################################################################################
// ###########################################################################################################
// ###########################################################################################################
void SWAP_FORWARD(MPS& psi, const std::vector<Index>& bin, int from, int to, double cutoff)
{
    ITensor SWAP,U,S,V;
    Index iSWAP;
    auto SWITCH_L = findIndex(psi(from),"left" );
    if ( SWITCH_L ) iSWAP=findIndex(psi(from),"left" );

    auto SWITCH_R = findIndex(psi(from),"right" );
    if ( SWITCH_R ) iSWAP=findIndex(psi(from),"right" );
         
    if ( SWITCH_L || SWITCH_R )
    {    
     for(int k=from;k<to;k++)
     {    
      SWAP = psi(k)*psi(k+1);
      U=ITensor(iSWAP,commonIndex(psi(k+1),psi(k+2)));
      svd(SWAP,U,S,V,{"Cutoff=",cutoff});
      psi.setA(k,V); 
      psi.setA(k+1,U*S);
     }
    }
    else { printf("ERROR!! \n"); };
}    
// ###########################################################################################################
// ###########################################################################################################
// ###########################################################################################################

// ###########################################################################################################
// ###########################################################################################################
// ###########################################################################################################
void SWAP_BACKWARD(MPS& psi, const std::vector<Index>& bin, int from, int to, double cutoff)
{
    ITensor SWAP,U,S,V;
    Index iSWAP;
    auto SWITCH_L = findIndex(psi(from),"left" );
    if ( SWITCH_L ) iSWAP=findIndex(psi(from),"left" );

    auto SWITCH_R = findIndex(psi(from),"right" );
    if ( SWITCH_R ) iSWAP=findIndex(psi(from),"right" );    
    
    for(int k=from;k>to;k--)
    {    
    SWAP = psi(k)*psi(k-1);
    if (order(psi(k-1))==2) U=ITensor(iSWAP); 
    else U=ITensor(iSWAP,commonIndex(psi(k-1),psi(k-2)));
    svd(SWAP,U,S,V,{"Cutoff=",cutoff});
    
    if (k-1==to) { psi.setA(k-1,U  ); psi.setA(k,V*S);}
    else         { psi.setA(k-1,U*S); psi.setA(k,V);  }
    }
}    
// ###########################################################################################################
// ###########################################################################################################
// ###########################################################################################################

// ###########################################################################################################
// ###########################################################################################################
// ###########################################################################################################
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
// ###########################################################################################################
// ###########################################################################################################
// ###########################################################################################################

int main(int argc, char* argv[])
{ 
if(argc != 2) { printfln("Usage: %s inputfile",argv[0]); return 0; } 
auto input = InputGroup(argv[1],"input");
// size of timestep
Real dt = input.getReal("time_step");
// maximum number of time steps
int t_end = input.getInt("time_end");
Real Gamma_l  = input.getReal("Gamma_l"); //Gamma_l = Gamma_l*sqrt(dt);
Real Gamma_ml = input.getReal("Gamma_ml"); //Gamma_m = Gamma_m*sqrt(dt);
Real Gamma_mr = input.getReal("Gamma_mr"); //Gamma_m = Gamma_m*sqrt(dt);
Real Gamma_r  = input.getReal("Gamma_r"); //Gamma_r = Gammm_r*sqrt(dt);
//dimension of local Hilbert space of each bin
int Nbin = input.getInt("Nbin",4); 
//cutoff of schmidt values
Real cutoff = input.getReal("svdcutoff");
//maximal number of schmidtvalues 
int maxm = input.getInt("maxnumberofSV");    
int fb  = input.getInt("fb");
int Nfb = 2*(fb+fb+fb); // size of array is twice because of left and right moving photons and one feedback loop 
int ttotal = 2*t_end+Nfb; printf("Total Number of time bins: %i\n",ttotal);
// ATTENTION at the moment I do not take phase in between the states into account
double Init[17];
Init[1]  = input.getReal("init_0000"); Init[2]  = input.getReal("init_0001"); 
Init[3]  = input.getReal("init_0010"); Init[4]  = input.getReal("init_0011"); 
Init[5]  = input.getReal("init_0100"); Init[6]  = input.getReal("init_0101"); 
Init[7]  = input.getReal("init_0110"); Init[8]  = input.getReal("init_0111");
Init[9]  = input.getReal("init_1000"); Init[10] = input.getReal("init_1001"); 
Init[11] = input.getReal("init_1010"); Init[12] = input.getReal("init_1011"); 
Init[13] = input.getReal("init_1100"); Init[14] = input.getReal("init_1101"); 
Init[15] = input.getReal("init_1110"); Init[16] = input.getReal("init_1111");
Real phi    = input.getReal("phi");
int  INDIVIDUAL_PHASE_CHOICE = input.getInt("INDIVIDUAL_PHASE_CHOICE");
Real phi_l  = input.getReal("phi_l");      
Real phi_ml = input.getReal("phi_ml");    
Real phi_mr = input.getReal("phi_mr");    
Real phi_r  = input.getReal("phi_r"); 
if ( (fb>0) && (INDIVIDUAL_PHASE_CHOICE == 0)) { phi_l=phi; phi_r=phi; phi_ml = phi*2./3.; phi_mr = phi*1./3.; }
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
double pop4=0.;
// ATOM 1  |9><9|  + |10><10| + |11><11| + |13><13| + |15><15| + |14><14| + |12><12| + |16><16|          
pop1 = Init[9 ]*Init[9 ]+Init[10]*Init[10]+Init[11]*Init[11]+Init[15]*Init[15]
      +Init[13]*Init[13]+Init[14]*Init[14]+Init[12]*Init[12]+Init[16]*Init[16] ; // ATTENTION Init_xyz is real
// ATOM 2 |5><5|   + |6><6|   + |7><7| + |13><13| + |15><15| + |14><14| + |8><8| + |16><16| 
pop2 = Init[5 ]*Init[5 ]+Init[6 ]*Init[6 ]+Init[7]*Init[7]+Init[13]*Init[13]
      +Init[15]*Init[15]+Init[14]*Init[14]+Init[8]*Init[8]+Init[16]*Init[16]   ; // ATTENTION Init_xyz is real
// ATOM § |3><3|   + |4><4|   + |7><7| + |11><11| + |15><15| + |12><12| + |8><8| + |16><16| 
pop3 = Init[3 ]*Init[3 ]+Init[4 ]*Init[4 ]+Init[7 ]*Init[7 ]+Init[11]*Init[11]
      +Init[15]*Init[15]+Init[12]*Init[12]+Init[8 ]*Init[8 ]+Init[16]*Init[16]; // ATTENTION Init_xyz is real
// ATOM 4 |2><2|   + |4><4|   + |6><6| + |10><10| + |14><14| + |12><12| + |8><8| + |16><16|    
pop4 = Init[2 ]*Init[2 ]+Init[4 ]*Init[4 ]+Init[6 ]*Init[6 ]+Init[10]*Init[10]
      +Init[14]*Init[14]+Init[12]*Init[12]+Init[8 ]*Init[8 ]+Init[16]*Init[16]; // ATTENTION Init_xyz is real

FILE *file;
FILE *f_prob;

char FILE_NAME[2048+1024];
snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_N4_OCCUPATIONS_fb%i_N%i_sum_pop_%.2f_p_l_%.2f_p_ml_%.2f_p_mr_%.2f_p_r_%.2f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,fb,t_end,pop4+pop3+pop2+pop1,phi_l,phi_ml,phi_mr,phi_r,dt,cutoff*1000000.);
file = fopen(FILE_NAME,"w");

snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_N4_PROBABILITY_fb%i_N%i_sum_pop_%.2f_p_l_%.2f_p_ml_%.2f_p_mr_%.2f_p_r_%.2f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,fb,t_end,pop4+pop3+pop2+pop1,phi_l,phi_ml,phi_mr,phi_r,dt,cutoff*1000000.);
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
int tls=0; bin.at(tls) = Index(16,"sys"); 
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

int one_l   = sys_at-2*(fb+fb+fb)   ; // Nfb=2*(r_fb+l_fb);
int two_l   = sys_at-2*(fb+fb   )   ;
int three_l = sys_at-2*(fb      )   ;
int four_l  = Nfb+1; 
int one_r   = Nfb+2;
int two_r   = sys_at-2*(fb      )+1 ;    
int three_r = sys_at-2*(fb+fb   )+1 ;    
int four_r  = sys_at-2*(fb+fb+fb)+1 ;

ITensor U_evo, temp; 
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

Real Prob[22]; 
Prob[20]=0.;
for(int p=1;p<=16;p++) {Prob[p]=state_prop(psi(sys_at),p); Prob[20] += Prob[p];}
               
ITensor U,S,V,W,SWAP;
Index iFBl,iFBr,iMRl,iMRr,iMLl,iMLr,iCBl,iCBr; // index fuer feedback, system, und current bin 
Index iLinkp; 
Index b; 
ITensor BdgB;

if (Nfb>3) // fb>1 between the emitter
{

MPO_SETUP(U_evo,bin,tls,one_l,two_l,three_l,four_l,one_r,two_r,three_r,four_r,Nbin,Gamma_l,Gamma_ml,Gamma_mr,Gamma_r,phi_l,phi_ml,phi_mr,phi_r,dt);

for(int m=0;m<t_end;m++)
{   
   one_l   = sys_at-2*(fb+fb+fb)   ; // Nfb=2*(r_fb+l_fb); atom 1 ... atom 2 ... atom 3 ... atom 4 
   four_r  = sys_at-2*(fb+fb+fb)+1 ;
   two_l   = sys_at-2*(fb+fb   )   ;
   three_r = sys_at-2*(fb+fb   )+1 ;    
   three_l = sys_at-2*(fb      )   ;
   two_r   = sys_at-2*(fb      )+1 ;    
   four_l  = sys_at                ; 
   one_r   = sys_at+1              ;

    // --- Status and File output ----------------------------------------------------------------------------------------
    sum_l +=l_out; sum_r += r_out;
    printf("Step %i of %i: norm=%.10f  -- ",m,t_end,mps_norm);
    printf("pop_l=%.10f -- pop_ml=%.10f -- pop_mr=%.10f -- pop_r=%.10f -- PROB=%.5f",pop1,pop2,pop3,pop4,Prob[20]); 
    //printf("l_out=%.10f -- r_out=%.10f -- Sum_Out=%.2f ",l_out,r_out,sum_l+sum_r); 
    printf("\n");
    fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \t %.10f \n",m*dt,pop1,pop2,pop3,pop4);
    fflush(file); 
    fprintf(f_prob,"%.10f \t",m*dt); 
    fprintf(f_prob,"%.10f \t",Prob[1]);
    fprintf(f_prob,"%.10f \t",Prob[2]);fprintf(f_prob,"%.10f \t",Prob[3]);fprintf(f_prob,"%.10f \t",Prob[5]); 
    fprintf(f_prob,"%.10f \t",Prob[4]);fprintf(f_prob,"%.10f \t",Prob[6]);fprintf(f_prob,"%.10f \t",Prob[7]); 
    fprintf(f_prob,"%.10f \t",Prob[8]);
    fprintf(f_prob,"%.10f \t",mps_norm);
    fprintf(f_prob,"%.10f \t %.10f \t",l_out,r_out);
    fprintf(f_prob,"\n");
    fflush(f_prob);
    // --- End - Status and File output -----------------------------------------------------------------------------------
    //PrintMPS(psi,1,ttotal+1);
    // -------------- PREPARE MPS TO APPLY MPO ------------------------------------------
    SWAP_FORWARD(psi,bin,two_r,sys_at-1,cutoff);
    SWAP_FORWARD(psi,bin,three_l,sys_at-1,cutoff);
    SWAP_FORWARD(psi,bin,three_r,sys_at-1,cutoff);
    SWAP_FORWARD(psi,bin,two_l,sys_at-1,cutoff);
    SWAP_FORWARD(psi,bin,four_r    ,sys_at-1,cutoff);
    SWAP_FORWARD(psi,bin,one_l    ,sys_at-1,cutoff);
    // --------------------------------------------------------
    // --------------------------------------------------------
    // ----------- APPLY MPO START ----------------------------
    iMLr = findIndex(psi.A(sys_at-6),"right"); 
    iMRl = findIndex(psi.A(sys_at-5),"left");      
    iMRr = findIndex(psi.A(sys_at-4),"right"); 
    iMLl = findIndex(psi.A(sys_at-3),"left");      
    iFBr = findIndex(psi.A(sys_at-2),"right"); 
    iFBl = findIndex(psi.A(sys_at-1),"left");      
    iCBl = findIndex(psi.A(sys_at+1),"left");
    iCBr = findIndex(psi.A(sys_at+2),"right");    
    iLinkp = commonIndex(psi.A(sys_at-7),psi.A(sys_at-6));
    temp = noPrime(psi(sys_at-6)*psi(sys_at-5)*psi(sys_at-4)*psi(sys_at-3)*psi(sys_at-2)*psi(sys_at-1)*U_evo*psi(sys_at)*psi(sys_at+1)*psi(sys_at+2));  
    //Print(temp);
    
    U=ITensor(iLinkp,iMLr,iMRl,iMRr,iMLl,iFBl,iFBr,iCBl,iCBr); 
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    // tls in the next mps slot 
    psi.setA(sys_at+2,V); 
    Prob[20]=0;
    for(int p=1;p<=16;p++) {Prob[p]=state_prop(V*S,p); Prob[20] += Prob[p]; }    
    pop1 = Prob[9 ]*Prob[9 ]+Prob[10]*Prob[10]+Prob[11]*Prob[11]+Prob[15]*Prob[15]
          +Prob[13]*Prob[13]+Prob[14]*Prob[14]+Prob[12]*Prob[12]+Prob[16]*Prob[16] ; 
    pop2 = Prob[5 ]*Prob[5 ]+Prob[6 ]*Prob[6 ]+Prob[7 ]*Prob[7 ]+Prob[13]*Prob[13]
          +Prob[15]*Prob[15]+Prob[14]*Prob[14]+Prob[8 ]*Prob[8 ]+Prob[16]*Prob[16] ; 
    pop3 = Prob[3 ]*Prob[3 ]+Prob[4 ]*Prob[4 ]+Prob[7 ]*Prob[7 ]+Prob[11]*Prob[11]
          +Prob[15]*Prob[15]+Prob[12]*Prob[12]+Prob[8 ]*Prob[8 ]+Prob[16]*Prob[16] ; 
    pop4 = Prob[2 ]*Prob[2 ]+Prob[4 ]*Prob[4 ]+Prob[6 ]*Prob[6 ]+Prob[10]*Prob[10]
          +Prob[14]*Prob[14]+Prob[12]*Prob[12]+Prob[8 ]*Prob[8 ]+Prob[16]*Prob[16] ; 
    mps_norm = eltC( dag(V*S)*V*S).real(); 
    // now factorize step by step the reservoir tensor moving orthoCenter
    temp=U*S; U=ITensor(iLinkp,iMLr,iMRl,iMRr,iMLl,iFBl,iFBr,iCBl); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at+1,V); // iCBr 
    temp=U*S; U=ITensor(iLinkp,iMLr,iMRl,iMRr,iMLl,iFBl,iFBr); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at,V); // iCBl 
    temp=U*S; U=ITensor(iLinkp,iMLr,iMRl,iMRr,iMLl,iFBl); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at-1,V); // iFBr
    // ---------------------------------------
    b=findIndex(V,"right");
    BdgB = ITensor(b,prime(b)); 
    BdgB.set(b(2),prime(b(2)),1.); // ATTENTION bath bin maximal 1 photon!!! Nbin 2!!! 
    r_out = eltC( dag(V*S)*noPrime(BdgB*V*S) ).real();
    // ---------------------------------------
    temp=U*S; U=ITensor(iLinkp,iMLr,iMRl,iMRr,iMLl); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at-2,V); // iFBl 
    // ---------------------------------------
    b=findIndex(V,"left");
    BdgB = ITensor(b,prime(b)); 
    BdgB.set(b(2),prime(b(2)),1.); // ATTENTION bath bin maximal 1 photon!!! Nbin 2!!! 
    l_out = eltC( dag(V*S) * noPrime(BdgB*V*S) ).real();
    // ---------------------------------------
    temp=U*S; U=ITensor(iLinkp,iMLr,iMRl,iMRr); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at-3,V); // iMLl 
    temp=U*S; U=ITensor(iLinkp,iMLr,iMRl); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at-4,V); // iMRr 
    temp=U*S; U=ITensor(iLinkp,iMLr); svd(temp,U,S,V,{"Cutoff=",cutoff}); psi.setA(sys_at-5,V); // iMRl 
    psi.setA(sys_at-6,U*S); //iMLr    
    // ----------- SWAP BACK START -----------------------------
    SWAP_BACKWARD(psi,bin,sys_at-6,two_r-5,cutoff);
    SWAP_BACKWARD(psi,bin,sys_at-5,three_l-4,cutoff);
    SWAP_BACKWARD(psi,bin,sys_at-4,three_r-3,cutoff);
    SWAP_BACKWARD(psi,bin,sys_at-3,two_l-2,cutoff);
    SWAP_BACKWARD(psi,bin,sys_at-2,one_l      ,cutoff);
    SWAP_BACKWARD(psi,bin,sys_at-1,four_r      ,cutoff); // 
    // ----------- SWAP BACK END ------------------------------
    // --------------------------------------------------------
    // Update MPO, only if a further loop is possible
    if (m<t_end-1)
    {       
    U_evo = delta(bin[four_l     ],bin[four_l+2     ])*U_evo*delta(prime(bin[four_l     ]),prime(bin[four_l+2     ])) ;
    U_evo = delta(bin[one_r     ],bin[one_r+2     ])*U_evo*delta(prime(bin[one_r     ]),prime(bin[one_r+2     ])) ;
    U_evo = delta(bin[three_l],bin[three_l+2])*U_evo*delta(prime(bin[three_l]),prime(bin[three_l+2])) ;
    U_evo = delta(bin[three_r],bin[three_r+2])*U_evo*delta(prime(bin[three_r]),prime(bin[three_r+2])) ;
    U_evo = delta(bin[two_l],bin[two_l+2])*U_evo*delta(prime(bin[two_l]),prime(bin[two_l+2])) ;
    U_evo = delta(bin[two_r],bin[two_r+2])*U_evo*delta(prime(bin[two_r]),prime(bin[two_r+2])) ;
    U_evo = delta(bin[one_l    ],bin[one_l+2    ])*U_evo*delta(prime(bin[one_l    ]),prime(bin[one_l+2    ])) ;
    U_evo = delta(bin[four_r    ],bin[four_r+2    ])*U_evo*delta(prime(bin[four_r    ]),prime(bin[four_r+2    ])) ;
    }
    // increase the position of the system by 2!!
    sys_at++;sys_at++;
}   
}

// #################################################################################################
// ################################### END fb > 2 ##################################################
// #################################################################################################

// #################################################################################################
// ################################### fb == 2 #####################################################
// #################################################################################################
// ------ if fb==2 no swapping, due to the orthocenter not to be combined with fb>2 ----------------

// #################################################################################################
// ################################### END fb == 2 ##################################################
// #################################################################################################

// #################################################################################################
// ################################### fb == 2 #####################################################
// #################################################################################################
if (Nfb==0)
{
    
MPO_NO_FB_SETUP(U_evo,bin,tls,four_l,Nbin,Gamma_l,Gamma_ml,Gamma_mr,Gamma_r,dt,phi_l,phi_ml,phi_mr,phi_r);
 
for(int m=0;m<t_end;m++)
{   
    four_l = sys_at        ; 
    one_r = sys_at+1;
    // --- Status and File output ----------------------------------------------------------------------------------------
    sum_l +=l_out; sum_r += r_out;
    printf("Step %i of %i: norm=%.10f  -- ",m,t_end,mps_norm);
    printf("pop_l=%.6f -- pop_ml=%.6f -- pop_mr=%.6f -- pop_r=%.10f --",pop1,pop2,pop3,pop4); 
    printf("l_out=%.10f -- r_out=%.10f -- Sum_Out=%.2f ",l_out,r_out,sum_l+sum_r); 
    printf("\n");
    fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \n",m*dt,pop1,pop2,pop3,pop4,mps_norm);
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
    // -------------- PREPARE MPS TO APPLY MPO ------------------------------------------
    iCBl = findIndex(psi.A(sys_at+1),"left"); //Print(psi(sys_at+1));
    iCBr = findIndex(psi.A(sys_at+2),"right"); //Print(psi(sys_at+2));    
    iLinkp = commonIndex(psi.A(sys_at),psi.A(sys_at-1)); //PrintData(U_evo);
    temp = noPrime(U_evo*psi(sys_at)*psi(sys_at+1)*psi(sys_at+2)); //Print(temp);
    U=ITensor(iLinkp,iCBl,iCBr);     
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    // tls in the next mps slot 
    psi.setA(sys_at+2,V); 
    Prob[20]=0;
    for(int p=1;p<=16;p++) {Prob[p]=state_prop(V*S,p); Prob[20] += Prob[p]; }
    pop1 = Prob[9 ]*Prob[9 ]+Prob[10]*Prob[10]+Prob[11]*Prob[11]+Prob[15]*Prob[15]
          +Prob[13]*Prob[13]+Prob[14]*Prob[14]+Prob[12]*Prob[12]+Prob[16]*Prob[16] ; 
    pop2 = Prob[5 ]*Prob[5 ]+Prob[6 ]*Prob[6 ]+Prob[7 ]*Prob[7 ]+Prob[13]*Prob[13]
          +Prob[15]*Prob[15]+Prob[14]*Prob[14]+Prob[8 ]*Prob[8 ]+Prob[16]*Prob[16] ; 
    pop3 = Prob[3 ]*Prob[3 ]+Prob[4 ]*Prob[4 ]+Prob[7 ]*Prob[7 ]+Prob[11]*Prob[11]
          +Prob[15]*Prob[15]+Prob[12]*Prob[12]+Prob[8 ]*Prob[8 ]+Prob[16]*Prob[16] ; 
    pop4 = Prob[2 ]*Prob[2 ]+Prob[4 ]*Prob[4 ]+Prob[6 ]*Prob[6 ]+Prob[10]*Prob[10]
          +Prob[14]*Prob[14]+Prob[12]*Prob[12]+Prob[8 ]*Prob[8 ]+Prob[16]*Prob[16] ; 
    mps_norm = eltC( dag(V*S)*V*S).real(); 
    // now factorize step by step the reservoir tensor moving orthoCenter
    temp=U*S; U=ITensor(iLinkp,iCBl); svd(temp,U,S,V,{"Cutoff=",cutoff}); 
    psi.setA(sys_at+1,V*S); // iCBr 
    psi.setA(sys_at,U); // iCBl 
    
    // ---------------------------------------
    b=findIndex(V,"right");
    BdgB = ITensor(b,prime(b)); 
    BdgB.set(b(2),prime(b(2)),1.); // ATTENTION bath bin maximal 1 photon!!! Nbin 2!!! 
    r_out = eltC( dag(V*S)*noPrime(BdgB*V*S) ).real();
    // ---------------------------------------
    b=findIndex(U,"left");
    BdgB = ITensor(b,prime(b)); 
    BdgB.set(b(2),prime(b(2)),1.); // ATTENTION bath bin maximal 1 photon!!! Nbin 2!!! 
    l_out = eltC( dag(U*S) * noPrime(BdgB*U*S) ).real();
    
    temp = psi(sys_at+1)*psi(sys_at+2);
    iLinkp=commonIndex(psi(sys_at+1),psi(sys_at));  
    U = ITensor(iCBr,iLinkp);
    svd(temp,U,S,V,{"Cutoff=",cutoff});
    psi.setA(sys_at+1,U);       
    psi.setA(sys_at+2,V*S);      
    // Update MPO, only if a further loop is possible
    if (m<t_end-1)
    {       
    U_evo = delta(bin[four_l],bin[four_l+2])*U_evo*delta(prime(bin[four_l]),prime(bin[four_l+2])) ;
    U_evo = delta(bin[one_r],bin[one_r+2])*U_evo*delta(prime(bin[one_r]),prime(bin[one_r+2])) ;
    }
    // increase the position of the system by 2!!
    sys_at++;sys_at++;
} // end time loop  
}
// #################################################################################################
// ################################### END fb == 0 #################################################
// #################################################################################################

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