// Willett's Capillary Force Model
// Willett et al., Langmuir 2000,16
    double f1,f2,f3,f4;

    double tmpVar  = 2.3*log(max(1e-64,lBond));
    double tmpVar2 = tmpVar*tmpVar;
    double tmpVar3 = tmpVar*tmpVar2;
    
    f1 = (-0.44507 ) 
       + (-0.1119)    * tmpVar
       + (-0.0121010) * tmpVar2
       + (-0.000500)  * tmpVar3;
       
    f2 = (1.9222)    
       + (-0.0668)    * tmpVar
       + (-0.0013375) * tmpVar2;
    
    f3 = (1.268)     
       + (0.19800)    * tmpVar
       + (0.02232000) * tmpVar2
       + (0.0008585)  * tmpVar3;
       
    f4 = (-0.010703) 
       + (0.03345)    * tmpVar 
       + (0.00185740) * tmpVar2;

    if (delta > 0)
    {
      double tmpVar_1= 2.3 *log(0.5*delta);
      Fn_coh = - 2.0 * M_PI 
                     * surfaceTension 
                     * sqrt(radi *radj) 
                     * (
                        exp(   
                              f1
                             -(
                                  f2 * exp(  
                                             f3 * tmpVar_1 
                                            +f4 * tmpVar_1 * tmpVar_1
                                          )
                              )
                           )
                       );
    }
    else if (delta <= 0)
    {
      Fn_coh = - 2.0 * M_PI 
                     * surfaceTension 
                     * sqrt(radi *radj) 
                     * exp(f1);
    }
