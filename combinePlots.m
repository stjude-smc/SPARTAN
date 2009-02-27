function h1 = combinePlots()
%MAKEPLOTS  Combines plots from multiple datasets
% 
%   Creates population histogram (_hist.txt) and TD plot (_tdp.txt) files
%   from all data within a directory.  This is done once for every
%   directory given by the user.  A plot (a la makeplots.m) is then created.
%   NOTE: all directories must have same number of such files.
% 

% TODO: combined state hists

% Prompt user for locations of each dataset to be combined
dirNames  = cell(0,1);
pop_files = cell(0,0);
tdp_files = cell(0,0);


disp('Select directories, hit cancel when finished');
while 1,
    % Get next directory name from user
    datapath = uigetdir();
    if datapath==0, break; end  %user hit "cancel"
    dirNames{end+1} = datapath;
    
    % Find all pop. histograms within directory, add to <pop_files>
    d = dir( [datapath filesep '*_hist.txt'] );
    pop_files(end+1,:) = strcat( [datapath filesep], {d.name} );
    
    % Find all TD plots within directory, add to <tdp_files>
    d = dir( [datapath filesep '*_tdp.txt'] );
    tdp_files(end+1,:) = strcat( [datapath filesep], {d.name} );
    
    d = dir([datapath filesep '*combined*']);
    assert( size(d,1)==0, 'Must remove old combined plots before starting!' );
end

nDirs  = size(pop_files,1);
nFiles = size(pop_files,2);


% Generate combined pop. histograms
for i=1:nDirs,  %for each directory
        
    %tot_hist = [];
    %nt = 0;
    
    % Generate combined population histogram
    for j=1:nFiles,  %for each pophist file
    
        cplotdata = load( pop_files{i,j} );
        
        % Create combined histogram if this is first file
        if j==1,
            tot_hist = cplotdata;
            nt = size(tot_hist,2); %number of frames+1
            
        % Otherwise, add the results to this one
        else
            % Resize histograms to of different length (nFrames)
            if size(cplotdata,2) > nt,
                cplotdata = cplotdata(:,1:nt);
            end
            %if size(cplotdata,2) < nt .....
            
            assert( sum(tot_hist(1,2:end) == cplotdata(1,2:end)) & ...
                    sum(tot_hist(2:end,1) == cplotdata(2:end,1)), ...
                    'Mismatched pop histogram axes' );
            
            tot_hist(2:end,2:end) = tot_hist(2:end,2:end) + cplotdata(2:end,2:end);
        end
        
    end
    
    % Save the combined pop histogram to file
    dirname = getdirname(dirNames{i});
    filename = [dirNames{i} filesep dirname '_combined_hist.txt'];
    dlmwrite( filename, tot_hist, ' ' );
    
end



% Generate combined TD plots
% TODO: use sum rather than average
for i=1:nDirs,  %for each directory
        
    %tot_tdp = [];
    %tot_time
    %nt = 0;
    
    % Generate combined TD plot
    for j=1:nFiles,  %for each pophist file
    
        tdpdata  = load( tdp_files{i,j} );
        tdp_time = tdpdata(1,1);
        assert( tdp_time ~=0, 'Time not saved in TD plot!' );
        
        % Create combined histogram if this is first file
        if j==1,
            tot_tdp = tdpdata;
            nt = size(tot_tdp,2); %number of frames+1
            tot_time = tdp_time;
            
        % Otherwise, add the results to this one
        else
            % Resize histograms to of different length (nFrames)
            if size(tdpdata,2) > nt,
                tdpdata  = tdpdata(:,1:nt);
            end
            %if size(tdpdata,2) < nt .....
            
            assert( sum(tot_tdp(1,2:end) == tdpdata(1,2:end)) & ...
                    sum(tot_tdp(2:end,1) == tdpdata(2:end,1)), ...
                    'Mismatched TD plot axes' );
            
            tot_tdp(2:end,2:end) = tot_tdp(2:end,2:end) + ...
                tdpdata(2:end,2:end)*tdp_time;
            tot_time = tot_time+tdp_time;
        end
        
    end
    
    % Save the combined pop histogram to file
    tot_tdp(2:end,2:end) = tot_tdp(2:end,2:end)/tot_time;  %normalize
    tot_tdp(1,1) = tot_time;
    
    dirname = getdirname(dirNames{i});
    filename = [dirNames{i} filesep dirname '_combined.qub_tdp.txt'];
    dlmwrite( filename, tot_tdp, ' ' );
    
end


disp('Finished.');



function output = getdirname( path )
% Gets Nth token in R (seperated by delim)

while numel(path)>0,
    [output,path] = strtok(path, filesep);
end








