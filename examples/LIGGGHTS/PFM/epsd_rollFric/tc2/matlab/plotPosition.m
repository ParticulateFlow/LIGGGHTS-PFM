%% clear workspace
clear all
close all
clc

%% ### user input #########################################################

% old style
% matFileName = 'data_cdt_origData.mat';
% dumpFiles = '../post/dump_*';
% override = false;

% infos from LIGGGHTS-input file
dumpfreqForce = 1;

dirfile = '..';
filepattern = 'data.*.txt';

% column in the force matrix
col_X = 1;
col_Z = 3;
col_tqY = 8;

% ATTENTION: Fieldnames are hardcoded. If filepattern changes, you have to
% adapt the code.

% #########################################################################

%% read all data

listfile = dir(fullfile(dirfile,filepattern));
data = readLogData(dirfile,listfile);
nFiles = length(listfile);

for ii=1:length(data)
    data(ii).timesteps = (1:1:length(data(ii).values(:,1))).*dumpfreqForce.*data(ii).dt;
end


%% init data
% init figures
hFig(1) = figure; % position
hFig(2) = figure; % torque
hFig(3) = figure; % plotyy

% legend parameters
leg = cell(4,1); % 4 figures
iCnt = ones(4,1);

cmap = colormap(lines(nFiles+2));
linS = {'-';'--';'-.';':'};

for ii=1:nFiles
    
    fname = data(ii).name;
    timesteps = data(ii).timesteps;
    
    % plot position and velocity
    % if size(data(ii).values,2)>=col_Z
    xpos = data(ii).values(:,col_X);
    xvel = (data(ii).values(2:end,col_X)-data(ii).values(1:end-1,col_X))/(dumpfreqForce*data(ii).dt);
    
%     zpos = data(ii).values(:,col_Z);
%     zvel = (data(ii).values(2:end,col_Z)-data(ii).values(1:end-1,col_Z))/(dumpfreqForce*data(ii).dt);
    
    figure(hFig(1));
    subplot(2,1,1); hold on
    plot(timesteps,xpos,'Color',cmap(ii,:));
    xlim([0 0.8]);
    
    subplot(2,1,2); hold on
    plot(timesteps(1:end-1),xvel,'Color',cmap(ii,:));
    
    % legend for postion plot
    leg{1,iCnt(1)} = fname;
    iCnt(1) = iCnt(1)+1;
    % end
    
    % plot torque
    tqy = data(ii).values(:,col_tqY);
    
    figure(hFig(2)); hold on
    plot(timesteps,tqy,'Color',cmap(ii,:));
    % legend for postion plot
    leg{2,iCnt(2)} = fname;
    iCnt(2) = iCnt(2)+1;
    
    % yyplot
    figure(hFig(3)); hold on
    clf
    [AX,h1,h2]=plotyy(timesteps,xpos,timesteps,tqy);
    %xlim(AX(1),[0.1 0.8]); xlim(AX(2),[0.1 0.8]);
    legend('Rolling distance','Total torque','Location','SouthEast');
    set(get(AX(1),'Ylabel'),'String','ylabel1') 
    set(get(AX(2),'Ylabel'),'String','ylabel2') 
    xlabel('xlabel');
    
    saveas(hFig(3),['tc2_plot_',data(ii).rfstyle],'epsc');
end


%% read data (old)

%clear data
%
% if (override || exist(matFileName,'file') ~= 2)
%     data = getDumpData(dumpFiles);
%     save(matFileName,'data','-mat');
% else
%     load(matFileName)
% end
%
% % from input script:
% % dump order .. [timestep nAtoms] id type x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius
% for i=1:size(data,2)
%     plotData.(data{1,i,1}) = cell2mat(data(:,i,2));
% end
% plotData.time = plotData.timestep.*dt;
%
% %% plot data
%
% figure;
% plot(plotData.time,plotData.x,plotData.time,plotData.y,plotData.time,plotData.z);
% legend('x','y','z');
%
% figure;
% plot(plotData.time,plotData.omegax,plotData.time,plotData.omegay,plotData.time,plotData.omegaz);
% legend('wx','wy','wz');
%
% figure;
% plot(plotData.time,plotData.vx,plotData.time,plotData.vy,plotData.time,plotData.vz);
% legend('vx','vy','vz');
%
% figure;
% plot(plotData.time,plotData.fx,plotData.time,plotData.fy,plotData.time,plotData.fz);
% legend('fx','fy','fz');
%
% figure;
% plot(plotData.time,plotData.tqx,plotData.time,plotData.tqy,plotData.time,plotData.tqz);
% legend('mx','my','mz');
%
% figure;
% plot(plotData.time,plotData.x./max(plotData.x),plotData.time,plotData.vx./max(plotData.vx),plotData.time,plotData.omegay./max(plotData.omegay));
% legend('x','vx','wy');