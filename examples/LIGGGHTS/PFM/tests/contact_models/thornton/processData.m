#!/usr/bin/octave --silent

% clear workspace
clear all
close all
clc

dirfile = './post';
filepattern = 'displacement*.txt';

% list all data files

filelist = dir(fullfile(dirfile,filepattern));

fout1 = fopen('./post/force_displacement.txt','w'); % open new data file
fout2 = fopen('./post/force_time.txt','w'); % open new data file

nFiles = length(filelist);

for ii=1:nFiles

    % name
    iName = filelist(ii).name;
    %disp(['Processing ',iName,' ...']);
    % values
    fin = fopen(fullfile(dirfile,iName), 'r'); % open file for reading

    % get 4th line to check if there is any contact data
    line = fgetl(fin); line = fgetl(fin);
    line = fgetl(fin); line = fgetl(fin);
    A = sscanf(line,'%d');
    if A>0
        % get 10th line to extract contact data
        line = fgetl(fin); line = fgetl(fin);
        line = fgetl(fin); line = fgetl(fin);
        line = fgetl(fin); line = fgetl(fin);
        A = sscanf(line,'%f',[1 9]);
        % write tangential overlap (micrometer) and Force (N) to the new file
        fprintf(fout1,'%f %f\n',A(8)*1000000,-A(5));
        % write time (microseconds) and Force (N) to the new file
        fprintf(fout2,'%f %f %f %f\n',(ii-1)*2-100,A(3),A(3)*0.1,-A(5));
    end
    fclose(fin);
end

fclose(fout1);
fclose(fout2);

