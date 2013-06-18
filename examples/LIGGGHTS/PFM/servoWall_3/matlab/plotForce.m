%% clear workspace
clear all
close all
clc

%% ### user input #########################################################
dirfile = '..';
filepattern = 'force.cad*';

% column in the force matrix
col_fZ = 3;
col_Z = 9;

% #########################################################################


listfile = dir(fullfile(dirfile,filepattern));
nFiles = length(listfile);

data = readLogData(dirfile,listfile);

cmap = colormap(lines(nFiles+2));

for ii=1:nFiles
    
    time = 0:data(ii).dt:data(ii).dt*(size(data(ii).values,1)-1);
    area = (2*data(ii).rad*data(ii).dcyldp)^2*0.25*pi;
    
    figure(1); hold on
    plot(time,data(ii).values(:,col_fZ)./area,'Color',cmap(ii,:));
    title(listfile(ii).name,'Interpreter','none');
    xlabel('Timesteps');
    ylabel('f\_total\_ in N');

    disp(['Mean for ',num2str(col_fZ),' direction = ',num2str(mean(data(ii).values(:,col_fZ)/area))]);
    
    figure(2); hold on
    subplot(2,1,1);
    plot(time,data(ii).values(:,col_Z),'Color',cmap(ii,:));
    subplot(2,1,2);
    plot(time(1:end-1),data(ii).values(2:end,col_Z)-data(ii).values(1:end-1,col_Z),'Color',cmap(ii,:));
end