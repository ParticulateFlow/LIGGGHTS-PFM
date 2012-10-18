clear all
close all
%clc

dirfile = '.';
listfile = dir(fullfile(dirfile,'force_cad1_f*'));


nFiles = length(listfile);
nColumn = 3;
nRow = ceil(nFiles/nColumn);

for j=1:nFiles
    
    disp(['Prozessing ',listfile(j).name,' ...']);

    tic();
    force = importdata(fullfile(dirfile,listfile(j).name)); % matlab-version
    %force = dlmread(fullfile(dirfile,listfile(j).name),' ',1,0); % octave dlmread is really slow!!!
    toc();
    
    for i=3 %1:size(force.data,2)
        subplot(nRow,nColumn,j);
        plot(force.data(:,i)); % matlab-version
        %plot(force(:,i)); % octave-version
        title(listfile(j).name,'Interpreter','none');
        xlabel('Timesteps');
        ylabel('f\_total\_ in N');

        %disp(['Mean for ',num2str(i),' direction = ',num2str(mean(force.data(:,i)))]); % does not work in octave
    end

    clear force
end
