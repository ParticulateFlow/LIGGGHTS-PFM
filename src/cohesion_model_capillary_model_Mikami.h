// Mikami's Capillary Force Model
// Mikami et al., CES 1998, 50(16)
    const double Aparam = -1.1*pow((max(1e-64,lBond)),-0.53);
    const double Bparam = -0.0082*log(max(1e-64,lBond))+0.48;
    const double Cparam = 0.0018*log(max(1e-64,lBond))+0.078;

    if (delta > 0)
      Fn_coh = - M_PI*surfaceTension*sqrt(radi*radj)*(exp(Aparam*delta+Bparam)+Cparam);
    else if (delta <= 0)
      Fn_coh = - M_PI*surfaceTension*sqrt(radi*radj)*(exp(Bparam)+Cparam);

