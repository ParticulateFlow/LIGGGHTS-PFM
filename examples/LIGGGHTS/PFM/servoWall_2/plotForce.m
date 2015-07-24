#!/usr/bin/octave --silent

%% clear workspace
clear all
close all
clc

%% ### user input #########################################################

deltat = 5e-5;
dumpfreq = 1e-3/deltat;

dirfile = './post';
filepattern = 'force.cad*';

% column in the force matrix
col_fZ = 3;
col_tZ = 6;
col_X = 7;
col_Z = 9;

saveDir = './post';

% save images
savePlot = true;

% save plot data
savePlotData = true;

% get from file name:
%           fieldname after underscore (lower case)
%           Number of case (Upper case or number)
expr = '((?<=[_.])[a-z]*)|((?<=[_.][a-z]*)([A-Z0-9]*(\.[e\-0-9]+)?))';


% area size according in.shearCell
cylArea = 0.025 .^2 * pi;

% ATTENTION: Fieldnames are hardcoded. If filepattern changes, you have to
% adapt the code.

% #########################################################################

%% read all data

listfile = dir(fullfile(dirfile,filepattern));
nFiles = length(listfile);

for ii=1:nFiles
    % name
    iName = listfile(ii).name;
    disp(['Processing ',iName,' ...']);
    % values
    force = importdata(fullfile(dirfile,iName), ' ');

    % get data from file name
    [flags] = regexp(iName,expr,'match');

    data(ii).name = iName;
    data(ii).header = force.textdata;
    data(ii).values = force.data;
    data(ii).timesteps = (1:1:length(data(ii).values(:,1))).*dumpfreq.*deltat;
    for jj=1:length(flags)/2
        data(ii).(flags{jj*2-1}) = flags{jj*2};
    end

    data(ii).area = cylArea;

end

%% init data
% init figures
hFig(1) = figure; % position
hFig(2) = figure; % z-force
hFig(3) = figure; % x-force

% legend parameters
leg = cell(3,1); % 3 figures

cmap = colormap(jet(nFiles));

for ii=1:nFiles

    fname = data(ii).name;
    timesteps = data(ii).timesteps;

    % plot z-vel
    if size(data(ii).values,2)>=col_Z && ...
            strcmp(data(ii).cad,'1')

        pos = data(ii).values(:,col_Z);
        vel = (data(ii).values(2:end,col_Z)-data(ii).values(1:end-1,col_Z))/(dumpfreq*deltat);

        figure(hFig(1));
        subplot(2,1,1); hold on
        plot(timesteps,pos,'Color',cmap(ii,:));

        subplot(2,1,2); hold on
        plot(timesteps(1:end-1),vel,'Color',cmap(ii,:));

        if savePlotData
            dlmwrite(fullfile(saveDir,['pos',num2str(ii),'.txt']), [ timesteps' , pos ], ' ');
            dlmwrite(fullfile(saveDir,['vel',num2str(ii),'.txt']), [ timesteps(1:end-1)' , vel ], ' ');
        end

        % legend for postion plot
        leg{1} = [leg{1}; fname];
    end

    % plot z-force
    if (strcmp(data(ii).cad,'1'))
        sigmaZ = data(ii).values(:,col_fZ);%./data(ii).area;
        mean_fZ = mean(sigmaZ);

        figure(hFig(2)); hold on
        plot(timesteps,sigmaZ,'Color',cmap(ii,:));

        if savePlotData
            dlmwrite(fullfile(saveDir,['forceZ',num2str(ii),'.txt']), [ timesteps' , sigmaZ ], ' ');
        end

        % legend for position plot
        leg{2} = [leg{2}; fname];
    end

    % plot total z-torque
    if (strcmp(data(ii).cad,'2'))
        tauXZ = data(ii).values(:,col_tZ);

        figure(hFig(3)); hold on
        plot(timesteps,tauXZ,'Color',cmap(ii,:));

        if savePlotData
            dlmwrite(fullfile(saveDir,['torqueZ',num2str(ii),'.txt']), [ timesteps' , tauXZ ], ' ');
        end

        %legend for force plots
        leg{3} = [leg{3}; fname];
    end
    
end

%% settings for figures
% generate legends
figure(hFig(1)); subplot(2,1,1);
hLeg = legend(leg{1});
set(hLeg,'Interpreter','none');
title('Position and Velocity','Interpreter','none','FontSize',16);
ylabel('position in m','FontSize',16);

subplot(2,1,2);
xlabel('time in s');
ylabel('velocity in m/s');

figure(hFig(2));
hLeg = legend(leg{2});
set(hLeg,'Interpreter','none','Location','SouthEast');
title('Force in z-direction','Interpreter','none');
xlabel('time in s');
ylabel('F_{z} in N');


figure(hFig(3));
hLeg = legend(leg{3});
set(hLeg,'Interpreter','none');
title('Torque','Interpreter','none');
xlabel('timesteps');
ylabel('M_{z} in Nm');

%% save figures
if savePlot
    for ii=1:length(hFig)
        print(hFig(ii), fullfile(saveDir,['figure',num2str(ii),'.svg']), '-dsvg');
        %print(hFig(ii), fullfile(saveDir,['figure',num2str(ii),'.eps']), '-deps');
    end
end
