%% clear workspace
clear all
close all
clc

%% ### user input #########################################################

matFileName = 'data_epsd.mat';
dumpFiles = '../post/dump_*';
dt = 1e-6;

% #########################################################################

%% read data

if (exist(matFileName,'file') ~= 2)
    data = getDumpData(dumpFiles);
    save(matFileName,'data','-mat');
else
    load(matFileName)
end

% from input script: 
% dump order .. [timestep nAtoms] id type x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius
timesteps = cell2mat(data(:,1,2));
time = timesteps.*dt;
x = cell2mat(data(:,5,2));
y = cell2mat(data(:,6,2));
z = cell2mat(data(:,7,2));
wx = cell2mat(data(:,17,2));
wy = cell2mat(data(:,18,2));
wz = cell2mat(data(:,19,2));

figure;
plot(time,x,time,y,time,z);
legend('x','y','z');
figure;
plot(time,wx,time,wy,time,wz);