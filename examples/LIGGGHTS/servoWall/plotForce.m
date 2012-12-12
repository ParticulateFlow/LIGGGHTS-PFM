%% clear workspace
clear all
close all
clc

%% ### user input #########################################################

deltat = 5e-5;
dumpfreq = 1e-3/deltat;

dirfile = '.';
filepattern = 'force_cad*_f*';

% column in the force matrix
col_fX = 1;
col_fZ = 3;
col_X = 7;
col_Z = 9;

% get from file name: 
%           fieldname after underscore (lower case)
%           Number of case (Upper case or number)
expr = '((?<=_)[a-z]*)|((?<=_[a-z]*)[A-Z0-9]*)';

% cylinder area
cylArea = 0.025 .^2 * pi;

% ATTENTION: Fieldnames are hardcoded. If filepattern changes, you have to
% adapt the code.

% compare experimental data
exp_flag = false;
exp_area = 0.05*0.05*pi;
exp_dir = '/media/sdb1/Projekte/_intern/materialProperties/poorManShearCell/experiment/500mue_material';
exp_file = 'Versuch06_500mue_1374g';%'Versuch04_500mue_484g';

testFlag = false;

% #########################################################################

%% read all data

listfile = dir(fullfile(dirfile,filepattern));
nFiles = length(listfile);

for ii=1:nFiles
    % name
    iName = listfile(ii).name;
    disp(['Prozessing ',iName,' ...']);
    % values
    force = importdata(fullfile(dirfile,iName));
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

% legend parameters
leg = cell(4,1); % 4 figures

cmap = colormap(lines(nFiles+2));

for ii=1:nFiles
    
    fname = data(ii).name;
    timesteps = data(ii).timesteps;
    
    % plot z-vel
    if size(data(ii).values,2)>=col_Z && ...
            strcmp(data(ii).cad,'1')
        
        if strcmp(data(ii).cad,'2')
            pos = data(ii).values(:,col_X);
            vel = (data(ii).values(2:end,col_X)-data(ii).values(1:end-1,col_X))/(dumpfreq*deltat);    
        else
            pos = data(ii).values(:,col_Z);
            vel = (data(ii).values(2:end,col_Z)-data(ii).values(1:end-1,col_Z))/(dumpfreq*deltat);
        end
        
        figure(hFig(1));
        subplot(2,1,1); hold on
        plot(timesteps,pos,'Color',cmap(ii,:));
        
        subplot(2,1,2); hold on
        plot(timesteps(1:end-1),vel,'Color',cmap(ii,:));
        
        % legend for postion plot
        leg{1} = [leg{1}; fname];
    end
    
    % plot z-force
    if (strcmp(data(ii).cad,'1'))
        sigmaZ = data(ii).values(:,col_fZ)./data(ii).area;
        mean_fZ = mean(sigmaZ);
        
        figure(hFig(2)); hold on
        plot(timesteps,sigmaZ,'Color',cmap(ii,:));
        
        % legend for postion plot
        leg{2} = [leg{2}; fname];
        
        disp([fname,': Mean for z direction = ',num2str(mean_fZ)]);
    end
    
end

%% load experimental data
if (exp_flag)
    load(fullfile(exp_dir,exp_file), 'Fr', 't');
    exp_shear = Fr./exp_area;

    figure(hFig(3)); hold on
    plot(t,exp_shear,'Color',cmap(nFiles+1,:),'LineWidth',2);
    
    %leg{3} = [leg{3}; exp_file];
end

%% settings for figures
% generate legends
figure(hFig(1)); subplot(2,1,1);
hLeg = legend(leg{1});
set(hLeg,'Interpreter','none');
title('Position and Velocity','Interpreter','none','FontSize',16);
ylabel('position in m','FontSize',16);

subplot(2,1,2);
xlabel('time in s'); %xlabel('timesteps');
ylabel('velocity in m/s');

figure(hFig(2));
hLeg = legend(leg{2});
set(hLeg,'Interpreter','none','Location','SouthEast');
title('Prinzipal stress in z-direction','Interpreter','none');
xlabel('time in s'); %xlabel('timesteps');
ylabel('sigma_{z} in Pa');

