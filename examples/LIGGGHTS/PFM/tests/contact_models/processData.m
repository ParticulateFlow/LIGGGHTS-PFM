#!/usr/bin/octave --silent

% clear workspace
clear all
close all
clc

dirfile = './post';
filepattern = 'displacement*.txt';

% list all data files

filelist = dir(fullfile(dirfile,filepattern));

fout = fopen('./post/force_displacement.txt','w'); % open new data file

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
        A = sscanf(line,'%f',[1 6]);
        % write tangential overlap (micrometer) and Force (N) to the new file
        fprintf(fout,'%f %f\n',A(5)*1000000,-A(2));
    end
    fclose(fin);
end

fclose(fout);

