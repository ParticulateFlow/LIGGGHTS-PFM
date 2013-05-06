%% clear workspace
clear all
close all
clc

%% ### user input #########################################################

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

    xpos = data(ii).values(:,col_X);
     
    figure(hFig(1)); hold on
    plot(timesteps,xpos,'Color',cmap(ii,:));
    xlim([0 0.8]);
    % legend for postion plot
    leg{1,iCnt(1)} = fname;
    iCnt(1) = iCnt(1)+1;

    
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
    set(get(AX(1),'Ylabel'),'String','ylabel1') 
    set(get(AX(2),'Ylabel'),'String','ylabel2') 
    xlabel('xlabel');
    
    saveas(hFig(3),['tc2_plot3_',data(ii).rfstyle],'epsc');
end

%% settings for plots

figure(hFig(1));
ylim([0.1046 0.10485]);
xlim([0.2 0.6]);
legend('CDT','EPSD','Location','SouthEast');
xlabel('xlabel');
ylabel('ylabel');

saveas(hFig(1),['tc2_plot1'],'epsc');
