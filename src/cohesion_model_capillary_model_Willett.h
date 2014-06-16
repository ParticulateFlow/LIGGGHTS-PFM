// Willett's Capillary Force Model
// Willett et al., Langmuir 2000,16
    double f1,f2,f3,f4;

    f1 = (-0.44507 ) + (-0.1119) *2.3 *log(lBond) + (-0.0121010) *pow(2.3 *log(lBond),2) + (-0.000500) *pow(2.3 *log(lBond),3);
    f2 = (1.9222)    + (-0.0668) *2.3 *log(lBond) + (-0.0013375) *pow(2.3 *log(lBond),2);
    f3 = (1.268)     + (0.19800) *2.3 *log(lBond) + (0.02232000) *pow(2.3 *log(lBond),2) + (0.0008585) *pow(2.3 *log(lBond),3);
    f4 = (-0.010703) + (0.03345) *2.3 *log(lBond) + (0.00185740) *pow(2.3 *log(lBond),2);

    if (delta > 0)
      Fn_coh = - 2.0 *M_PI *surfaceTension *sqrt(radi *radj) *(exp(f1-(f2 *exp(f3 *2.3 *log(delta/2) + f4 *pow(2.3 *log(delta/2),2)))));
    else if (delta <= 0)
      Fn_coh = - 2.0 *M_PI *surfaceTension *sqrt(radi *radj) *(exp(f1));
