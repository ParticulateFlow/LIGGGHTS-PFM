%% clear workspace
clear all
close all
clc

%% ### user input #########################################################

deltat = 5e-5;
dumpfreq = 1e-3/deltat;

dirfile = '.';
filepattern = 'force.cad*';

% column in the force matrix
col_fX = 1;
col_fZ = 3;
col_tZ = 6;
col_X = 7;
col_Z = 9;

% get from file name: 
%           fieldname after underscore (lower case)
%           Number of case (Upper case or number)
%expr = '((?<=_)[a-z]*)|((?<=_[a-z]*)[A-Z0-9]*)';
expr = '((?<=[_.])[a-z]*)|((?<=[_.][a-z]*)([A-Z0-9]*(\.[e\-0-9]+)?))';


% area flag correspond to area size according in.shearCell
% rX with X = [A B C D]  ==> cylRadius = 0.05*[0.25 0.5 0.75 1}
% at the moment hardcoded switch; check code!!
cylArea = (66.7e-3)^2*pi - (27.3e-3)^2*pi;%(0.05*0.25)^2 *pi;%(0.05.*[0.25 0.5 0.75 1]).^2 * pi;

% ATTENTION: Fieldnames are hardcoded. If filepattern changes, you have to
% adapt the code.

% compare experimental data
exp_flag = false;
exp_area = 0.05*0.05*pi;
exp_dir = '500mue_material';
exp_file = 'Versuch06_500mue_1374g';%'Versuch04_500mue_484g';

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
    
    data(ii).r = 'A';
    
    switch data(ii).r
        case 'A'; data(ii).area = cylArea(1);
        case 'B'; data(ii).area = cylArea(2);
        case 'C'; data(ii).area = cylArea(3);
        case 'D'; data(ii).area = cylArea(4);
        otherwise; error('ERROR');
    end
end

%% init data
% init figures
hFig(1) = figure; % position
hFig(2) = figure; % z-force
hFig(3) = figure; % x-force

% legend parameters
leg = cell(3,1); % 3 figures

cmap = colormap(hsv(nFiles+1));

for ii=1:nFiles
    
    fname = data(ii).name;
    timesteps = data(ii).timesteps;
    
    % plot z-vel
    if size(data(ii).values,2)>=col_Z && ...
            (strcmp(data(ii).cad,'2') || strcmp(data(ii).cad,'3'))
        
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
    if (strcmp(data(ii).cad,'3'))
        sigmaZ = data(ii).values(:,col_fZ)./data(ii).area;
        mean_fZ = mean(sigmaZ);
        
        figure(hFig(2)); hold on
        plot(timesteps,sigmaZ,'Color',cmap(ii,:));
        
        % legend for postion plot
        leg{2} = [leg{2}; fname];
        
        disp([fname,': Mean for z direction = ',num2str(mean_fZ)]);
    end
    
    % plot total x-force
    tauXZ = data(ii).values(:,col_tZ);%./data(ii).area/((66.7e-3+27.3e-3)*0.5);
    
    figure(hFig(3)); hold on
    plot(timesteps,tauXZ,'Color',cmap(ii,:));
    
    %legend for force plots
    leg{3} = [leg{3}; fname];
    
end

%% load experimental data
if (exp_flag)
    load(fullfile(exp_dir,exp_file), 'Fr', 't');
    exp_shear = Fr./exp_area;

    figure(hFig(3)); hold on
    plot (t,exp_shear,'Color',cmap(nFiles+1,:));
    
    %leg{3} = [leg{3}; exp_file];
end

%% settings for figures
% generate legends
figure(hFig(1)); subplot(2,1,1);
hLeg = legend(leg{1});
set(hLeg,'Interpreter','none');
title('Position and Velocity','Interpreter','none');
ylabel('position in m');

subplot(2,1,2);
xlabel('timesteps');
ylabel('velocity in m/s');

figure(hFig(2));
hLeg = legend(leg{2});
set(hLeg,'Interpreter','none');
title('Normal stress in z-direction','Interpreter','none');
xlabel('timesteps');
ylabel('sigma_{z} in Pa');


figure(hFig(3));
hLeg = legend(leg{3});
set(hLeg,'Interpreter','none');
title('Shear stress in x-direction','Interpreter','none');
xlabel('timesteps');
ylabel('tau_{xz} in Pa');

%% sum up x-forces

flag_done = zeros(nFiles,1);
for ii=1:nFiles
    
end