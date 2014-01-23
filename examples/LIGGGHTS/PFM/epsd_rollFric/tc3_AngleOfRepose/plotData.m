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
col_ke = 2;
col_rke = 3;

% ATTENTION: Fieldnames are hardcoded. If filepattern changes, you have to
% adapt the code.

% #########################################################################

%% read all data

listfile = dir(fullfile(dirfile,filepattern));
data = readLogData(dirfile,listfile);
nFiles = length(listfile);

%% init data
% init figures
hFig(1) = figure; % energy

% legend parameters
leg = cell(4,1); % 4 figures
iCnt = ones(4,1);

cmap = colormap(jet(nFiles+2));
linS = {'-';'--';'-.';':'};

for ii=1:nFiles

    fname = data(ii).name;
    timesteps = data(ii).values(:,col_t);

    % plot energy
    ke = data(ii).values(:,col_ke);
    rke = data(ii).values(:,col_rke);

    keRef = 1;%8/15*(5e-3)^5*pi*1056 * 0.5 * (180/pi)^2; % ke = I * 0.5 * w^2

    figure(hFig(1)); hold on
    step = 1;
    plot(timesteps(1:step:end),ke(1:step:end)./keRef,'Color',cmap(ii,:));
    set(gca,'YScale','Log');

    dlmwrite(fullfile(dirfile,['tc3_plot_',num2str(ii),'.txt']), [ timesteps, ke, rke ], ' ');
end

%% legend
figure(hFig(1));
%hLeg = legend(leg{1,1:iCnt(1)-1});
%set(hLeg,'Interpreter','none');
legend('EPSD','EPSD2','CDT','Location','SouthWest');
xlabel('time [s]');
ylabel('Ekin [J]');

%saveas(hFig(1),['tc3_plot1'],'epsc');
print(hFig(1), fullfile(dirfile, 'tc3_plot1.svg'), '-dsvg');


