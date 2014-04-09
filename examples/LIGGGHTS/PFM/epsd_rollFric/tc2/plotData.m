#!/usr/bin/octave --silent

%% clear workspace
clear all
close all
clc

%% ### user input #########################################################

dirfile = './post';
filepattern = 'data.*.txt';

% column in the force matrix
col_t = 1;
col_X = 2;
col_Z = 4;
col_tqY = 9;

% ATTENTION: Fieldnames are hardcoded. If filepattern changes, you have to
% adapt the code.

% #########################################################################

%% read all data

listfile = dir(fullfile(dirfile,filepattern));
data = readLogData(dirfile,listfile);
nFiles = length(listfile);

%% init data
% init figures
hFig(1) = figure; % position
hFig(2) = figure; % torque
hFig(3) = figure; % plotyy

% legend parameters
leg = cell(4,1); % 4 figures
iCnt = ones(4,1);

cmap = colormap(jet(nFiles+2));
linS = {'-';'--';'-.';':'};

for ii=1:nFiles

    fname = data(ii).name;
    timesteps = data(ii).values(:,col_t);

    % plot position and velocity
    xpos = data(ii).values(:,col_X);

    figure(hFig(1)); hold on
    hPlot(ii,1)=plot(timesteps(end/3:end),xpos(end/3:end),'Color',cmap(ii,:));
    hPlot(ii,2)=line([0.2 0.8],[xpos(end) xpos(end)],'Color',cmap(ii+1,:),'LineStyle','-.');
    xlim([0 0.8]);
    % legend for postion plot
    leg{1,iCnt(1)} = fname;
    iCnt(1) = iCnt(1)+1;

    cline(1:length(xpos)) = xpos(end);
    dlmwrite(fullfile(dirfile,['tc2_plot1_',num2str(ii),'.txt']), [ timesteps(end/3:end), xpos(end/3:end), cline(end/3:end)' ], ' ');

    % plot torque
    tqy = data(ii).values(:,col_tqY);

    figure(hFig(2)); hold on
    plot(timesteps,tqy,'Color',cmap(ii,:));
    % legend for torque plot
    leg{2,iCnt(2)} = fname;
    iCnt(2) = iCnt(2)+1;

    % yyplot
    figure(hFig(3)); hold on
    clf
    [AX,h1,h2]=plotyy(timesteps,xpos,timesteps,tqy);
    %xlim(AX(1),[0.1 0.8]); xlim(AX(2),[0.1 0.8]);
    legend('Rolling distance','Total torque','Location','SouthEast');
    set(get(AX(1),'Ylabel'),'String','x-position [m]')
    set(get(AX(2),'Ylabel'),'String','torque [Nm]')
    xlabel('time [s]');

    %saveas(hFig(3),['tc2_plot3_',data(ii).rfstyle],'epsc');
    print(hFig(3), fullfile(dirfile, ['tc2_plot3_',num2str(ii),'.svg']), '-dsvg');

    dlmwrite(fullfile(dirfile,['tc2_plot3_',num2str(ii),'.txt']), [ timesteps, xpos, tqy ], ' ');
end

%% settings for plots

figure(hFig(1));
ylim([0.1045 0.10475]);
xlim([0.2 0.8]);
legend(hPlot(:,1),'EPSD','EPSD2','CDT','Location','SouthEast');
xlabel('time [s]');
ylabel('x-position [m]');

%saveas(hFig(1),['tc2_plot1'],'epsc');
print(hFig(1), fullfile(dirfile, 'tc2_plot1.svg'), '-dsvg');

