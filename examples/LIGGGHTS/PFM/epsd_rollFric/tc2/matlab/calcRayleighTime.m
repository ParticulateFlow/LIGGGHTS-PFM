%% clear workspace
clear all
close all
clc

%% ### user input #########################################################

r =     5e-3;
rho =   1056;
nu =    0.49;
Y =     40e6;%65e9;

% #########################################################################

%% calc time

G = 1/((2*(2+nu)*(1-nu))/(Y));

dt_r = pi*r*sqrt(rho/G)/(0.1631*nu+0.8766);

disp(['The Rayleigh time step is ',num2str(dt_r),'. 10% of it are ',num2str(0.1*dt_r),'.']);