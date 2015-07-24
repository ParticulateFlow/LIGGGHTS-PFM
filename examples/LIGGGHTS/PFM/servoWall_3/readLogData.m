#!/usr/bin/octave --silent

function [ data ] = readLogData( dirfile,filelist )
%READLOGDATA Summary of this function goes here
%   Detailed explanation goes here

% get from file name:
%           fieldname after underscore (lower case)
%           Number of case (Upper case or number)
expr = '((?<=[_.])[a-z]*)|((?<=[_.][a-z]*)([A-Z0-9]*(\.[0-9]+)?(e\-[0-9]+)?))';

nFiles = length(filelist);

for ii=1:nFiles
    % name
    iName = filelist(ii).name;
    disp(['Processing ',iName,' ...']);
    % values
    rdata = importdata(fullfile(dirfile,iName), ' ');

    % get data from file name
    [flags] = regexp(iName,expr,'match');

    data(ii).name = iName;
    data(ii).header = rdata.textdata;
    data(ii).values = rdata.data;
    for jj=1:(length(flags)-1)/2
        if (isempty(regexp(flags{jj*2},'[A-Z]+'))) % it is a number
            data(ii).(flags{jj*2-1}) = str2num(flags{jj*2});
        else
            data(ii).(flags{jj*2-1}) = flags{jj*2};
        end
    end

end

end

