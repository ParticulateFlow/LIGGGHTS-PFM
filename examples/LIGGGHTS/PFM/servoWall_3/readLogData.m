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
    % rdata = importdata(fullfile(dirfile,iName)); % matlab-version
	% octave version
    rdata.data = [];
    rdata.textdata = [];

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
    rdata.textdata = file_content_rows(1:1)';

    data_columns = length(strsplit(file_content_rows{2}, " ", true)); % expecting all rows to have same # of columns

    rdata.data = NaN(length(file_content_rows) - 1, data_columns);
    for i=2:length(file_content_rows)
        % Only use the row if it contains anything other than white-space characters.
        if (any (file_content_rows{i} != ' '))
            [row_data] = strsplit(file_content_rows{i}, " ", true);

            for j=1:length(row_data)
                % Try to convert the column to a number, if it works put it in output.data, otherwise in output.textdata
                data_numeric = str2double (row_data{j});
                rdata.data(i-1, j) = data_numeric;
            endfor
        endif
    endfor
	% end octave version

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

