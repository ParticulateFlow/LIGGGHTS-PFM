#!/usr/bin/octave --silent

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
    disp(['Prozessing ',iName,' ...']);
    % values
    % force = importdata(fullfile(dirfile,iName)); % matlab version
	% octave version
    force.data = [];
    force.textdata = [];

    fid = fopen (fullfile(dirfile,iName), "r");
    if (fid < 0)
        error ("cannot open file");
    endif
    unwind_protect
    file_content = (fread (fid, '*char')).';
    unwind_protect_cleanup
    fclose (fid);
    end_unwind_protect

    [file_content_rows] = strsplit(file_content, "\n", true);
    force.textdata = file_content_rows(1:1)';

    data_columns = length(strsplit(file_content_rows{2}, " ", true)); % expecting all rows to have same # of columns

    force.data = NaN(length(file_content_rows) - 1, data_columns);
    for i=2:length(file_content_rows)
        % Only use the row if it contains anything other than white-space characters.
        if (any (file_content_rows{i} != ' '))
            [row_data] = strsplit(file_content_rows{i}, " ", true);

            for j=1:length(row_data)
                % Try to convert the column to a number, if it works put it in output.data, otherwise in output.textdata
                data_numeric = str2double (row_data{j});
                force.data(i-1, j) = data_numeric;
            endfor
        endif
    endfor
    % end octave version

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
%hFig(3) = figure; % x-force

% legend parameters
leg = cell(2,1); % 2 figures

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
        sigmaZ = data(ii).values(:,col_fZ)./data(ii).area;
        mean_fZ = mean(sigmaZ);

        figure(hFig(2)); hold on
        plot(timesteps,sigmaZ,'Color',cmap(ii,:));

        if savePlotData
            dlmwrite(fullfile(saveDir,['sigmaZ',num2str(ii),'.txt']), [ timesteps' , sigmaZ ], ' ');
        end

        % legend for position plot
        leg{2} = [leg{2}; fname];
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
title('Principal stress in z-direction','Interpreter','none');
xlabel('time in s');
ylabel('sigma_{z} in Pa');

%% save figures
if savePlot
    for ii=1:length(hFig)
        print(hFig(ii), fullfile(saveDir,['figure',num2str(ii),'.svg']), '-dsvg');
        %print(hFig(ii), fullfile(saveDir,['figure',num2str(ii),'.eps']), '-deps');
    end
end
