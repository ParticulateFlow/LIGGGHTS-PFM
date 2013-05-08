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
col_ke = 1;

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
hFig(1) = figure; % energy

% legend parameters
leg = cell(4,1); % 4 figures
iCnt = ones(4,1);

cmap = colormap(lines(nFiles+2));
linS = {'-';'--';'-.';':'};

for ii=1:nFiles
    
    fname = data(ii).name;
    timesteps = data(ii).timesteps;
    
    % plot energy
    ke = data(ii).values(:,col_ke);
    
    keRef = 1;%8/15*(5e-3)^5*pi*1056 * 0.5 * (180/pi)^2; % ke = I * 0.5 * w^2

    figure(hFig(1)); hold on
    step = 1;
    axe = plot(timesteps(1:step:end),ke(1:step:end)./keRef,'Color',cmap(ii,:));
    set(gca,'YScale','Log');

    % legend for postion plot
    leg{1,iCnt(1)} = fname;
    iCnt(1) = iCnt(1)+1;
    
end

%% legend
figure(hFig(1));
%hLeg = legend(leg{1,1:iCnt(1)-1});
%set(hLeg,'Interpreter','none');
legend('CDT','EPSD','Location','SouthWest');
xlabel('xlabel');
ylabel('ylabel');

saveas(hFig(1),['tc3_plot1'],'epsc');



