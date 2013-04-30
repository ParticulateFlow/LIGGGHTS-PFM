function [ data ] = getDumpData( fPath )
%GETDUMPDATA Read LIGGGHTS dump-files and returns data as cell array.
%   DATA = SUM(fPATH) is the cell array containing the values from the
%   dump-files at fPATH.
%
%   Returns the cell array DATA
%   DATA{idx_timestep, idx_property, 1/2}
%       1 ... name of property
%       2 ... value[s]
%
%   Examples:
%   DATA = getDumpData('./post/dump*.liggghts')
%
%   ATTENTION: No error checks are implemented.
%              It should work for both, Matlab and Octave.
%
%   Author: Andreas Aigner, Aug. 2012

disp('Read in data ...');

% list of all files
listFile = dir(fPath);
pathstr = fileparts(fPath);

for t=1:length(listFile)
    
    disp(listFile(t).name);
    fid = fopen(fullfile(pathstr,listFile(t).name));
    
    idx_prop = 1;
    
    tline = fgetl(fid);
    while not(feof(fid)) && ischar(tline)
        if isempty(tline)
            tline = fgetl(fid);
            continue;
        else
            if strcmp('ITEM: TIMESTEP',tline)
                % disp('TIMESTEP');
                dt = fscanf(fid,'%d \n');
                data{t,idx_prop,1} = 'timestep';
                data{t,idx_prop,2} = dt;
                idx_prop = idx_prop+1;
            elseif strcmp('ITEM: NUMBER OF ATOMS',tline)
                % disp('NUMBER OF ATOMS');
                n = fscanf(fid,'%d \n');
                data{t,idx_prop,1} = 'nAtoms';
                data{t,idx_prop,2} = n;
                idx_prop = idx_prop+1;
            elseif strcmp('ITEM: BOX BOUNDS',tline(1:16))
                % disp('BOX');
                box = fscanf(fid,'%f %f\n',[2 3]); % we ignore these values
            elseif strncmp('ITEM: ATOMS',tline,11)
                % disp('REST');
                % get names
                % id type type x y z ix iy iz vx vy vz fx fy fz vx vy vz p p rho
                tags = regexp(tline,'[a-z]+','match'); % return only lowercased properties
                values = fscanf(fid,'%f\n',[length(tags) n]); % get all values
                %values = fscanf(fid,'%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n',[21 n]);
                
                for i=0:length(tags)-1
                    data{t,idx_prop+i,1} = tags(i+1);
                    data{t,idx_prop+i,2} = values(i+1,:);
                end
                
                idx_prop = idx_prop + length(tags);
                %rho = values(21,:);
            end
        end
        
        tline = fgetl(fid);
    end
    
    clear dt n box tags values
    
    fclose(fid);
    
end

clear i t fid tline
disp('Finished.');

%% sort time data
disp('Sort time data ...');

% first sort algorithm (kind of a bubble sort)

% nTime = size(data,1);
% i = 1;
% change_flag = false;
% while i <= nTime
%     idx = i;
%     % disp(data{idx,1,2});
%     if (idx+1 <= nTime && data{idx,1,2} > data{idx+1,1,2})
%         %disp([num2str(data{idx,1,2}),' is bigger than ',num2str(data{idx+1,1,2})]);
%         swapIdx = idx+1;
%         change_flag = true;
%         while (swapIdx+1 <= nTime  && data{idx,1,2} > data{swapIdx+1,1,2})
%             %disp([num2str(data{idx,1,2}),' is bigger than ',num2str(data{swapIdx+1,1,2})]);
%             swapIdx = swapIdx+1;
%         end
%         
%         % shift data
%         tmp = data(idx,:,:);
%         tmpIdx = idx:swapIdx;
%         for j=2:length(tmpIdx)
%             data(tmpIdx(j-1),:,:) = data(tmpIdx(j),:,:);
%         end
%         data(swapIdx,:,:) = tmp;
%     end
%     
%     if (idx-1 ~= 0 && data{idx,1,2} < data{idx-1,1,2})
%         %disp([num2str(data{idx,1,2}),' is smaller than ',num2str(data{idx-1,1,2})]);
%         swapIdx = idx-1;
%         change_flag = true;
%         while (swapIdx-1 ~= 0 && data{idx,1,2} < data{swapIdx-1,1,2})
%             %disp([num2str(data{idx,1,2}),' is smaller than ',num2str(data{swapIdx-1,1,2})]);
%             swapIdx = swapIdx-1;
%         end
%         
%         % shift data
%         tmp = data(idx,:,:);
%         tmpIdx = idx:-1:swapIdx;
%         for j=2:length(tmpIdx)
%             data(tmpIdx(j-1),:,:) = data(tmpIdx(j),:,:);
%         end
%         data(swapIdx,:,:) = tmp;
%     end
% 
%     
%     if change_flag
%         change_flag = false;
%     else
%         i = i+1;
%     end
% end
% 
% clear change_flag tmp nTime i idx

% quick sort
set(0,'RecursionLimit',max(size(data,1),10000)); % Caution: This may crush your matlab!
quicksort(1,size(data,1));

disp('Finished.');


    % quick sort functions
    function [] = quicksort(left,right)
        if (left < right)
            divider = divide(left,right);
            quicksort(left,divider-1);
            quicksort(divider+1,right);
        end
    end
    function [posPivot] = divide(left,right)
        k = left;
        j = right-1;
        pivot = data{right,1,2};
        
        while (k<j)
            while (data{k,1,2} <= pivot && k < right)
                k = k+1;
            end
            
            while (data{j,1,2} >= pivot && j > left)
                j = j-1;
            end
            
            if (k<j)
                tmp = data(k,:,:);
                data(k,:,:) = data(j,:,:);
                data(j,:,:) = tmp;
            end
        end
        
        if (data{k,1,2} > pivot)
            tmp = data(k,:,:);
            data(k,:,:) = data(right,:,:);
            data(right,:,:) = tmp;
        end
        
        posPivot = k;
        
    end


end
