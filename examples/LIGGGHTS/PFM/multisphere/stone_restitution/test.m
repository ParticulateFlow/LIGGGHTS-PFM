#   mass = 3.112935e-03, radius of bounding sphere = 1.000000e-02, radius of equivalent sphere = 6.673912e-03
#   center of mass = 0.000000e+00, 0.000000e+00, 0.000000e+00
#   Principal moments of inertia: 3.377034e-08, 9.602173e-08, 9.595894e-08


v0 = 1.5
omegaY = 141;

m = 0.00311;
r = 0.005;
l = 0.02;

Ja = 0.5*m*r*r

Jb = 0.25*m*r*r + 1/12*m*l*l

Ekin = 0.5*m*v0*v0

Erot = 0.5*Jb*omegaY*omegaY




