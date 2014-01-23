#!/usr/bin/octave --silent

%% clear workspace
clear all
close all
clc

%% ### user input #########################################################
dumpfreq = 100;

dirfile = './post';
filepattern = 'force.cad*';

% column in the force matrix
col_fZ = 3;
col_Z = 9;

% #########################################################################


listfile = dir(fullfile(dirfile,filepattern));
nFiles = length(listfile);

data = readLogData(dirfile,listfile);

% init figures
hFig(1) = figure;
hFig(2) = figure;

cmap = colormap(jet(nFiles+2));

for ii=1:nFiles

    time = 0:data(ii).dt:data(ii).dt*(size(data(ii).values,1)-1);
    area = (2*data(ii).rad*data(ii).dcyldp)^2*0.25*pi;
    sigmaZ = data(ii).values(:,col_fZ)./area;
    pos = data(ii).values(:,col_Z);
    vel = (data(ii).values(2:end,col_Z)-data(ii).values(1:end-1,col_Z))/(dumpfreq*data(ii).dt);


    figure(hFig(1)); hold on
    plot(time,sigmaZ,'Color',cmap(ii,:));
    title(listfile(ii).name,'Interpreter','none');
    xlabel('Timesteps');
    ylabel('stress\_total\_ in Pa');

    print(hFig(1), fullfile(dirfile, 'figure1.svg'), '-dsvg');
    dlmwrite(fullfile(dirfile,['sigmaZ',num2str(ii),'.txt']), [ time' , sigmaZ ], ' ');

    disp(['Mean for ',num2str(col_fZ),' direction = ',num2str(mean(data(ii).values(:,col_fZ)/area))]);

    figure(hFig(2)); hold on
    subplot(2,1,1);
    plot(time,pos,'Color',cmap(ii,:));
    subplot(2,1,2);
    plot(time(1:end-1),data(ii).values(2:end,col_Z)-data(ii).values(1:end-1,col_Z),'Color',cmap(ii,:));

    dlmwrite(fullfile(dirfile,['pos',num2str(ii),'.txt']), [ time' , pos ], ' ');
    dlmwrite(fullfile(dirfile,['vel',num2str(ii),'.txt']), [ time(1:end-1)' , vel ], ' ');
end
