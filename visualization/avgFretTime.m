function varargout = avgFretTime( files )
%avgFretTime  Average FRET trajectory across all traces with bleaching removed.
%
%   avgFretTime(FILES) loads each .traces file in the cell array FILES, 
%   takes the median the FRET values from all traces in each frame, and plots
%   an average trajectory for each.
%
%   avgFretTime() will prompt the user for files.
%
%   OUT = avgFretTime(...) will save the timecourses in the matrix OUT with the
%   first column having the time axis in seconds. Nothing will be plotted.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


% Settings
truncateLen = 300;  %frames to calculate over
constants.min_fret = 0.175;  % minimum fret value, below which we assume there is no FRET.

% Get list of files if not specified.
if nargin<1, files = getFiles(); end
if isempty(files), return; end
if ~iscell(files), files = {files}; end
nFiles = numel(files);

% Generate plot titles and base filenames if none given.
[~,f] = cellfun(@fileparts, files, 'UniformOutput',false);
titles = trimtitles(f);



%% Calculate average FRET trajectories
output = zeros(truncateLen,1+nFiles);

for i=1:numel(files),
    % Load FRET data and truncate to target length
    data = loadTraces( files{i} );
    data.nFrames = truncateLen;
    
    if i==1,
        output(:,1) = data.time;  %convert to seconds
    else
        if ~all(data.time==output(:,1)'),
            warning('Time axes do not match!');
        end
    end
    
    % For each trace, average the FRET values at each frame to create an
    % average FRET trace. Exclude regions where the dyes are dark.
    for j=1:truncateLen,
        nonzero = data.fret(:,j) >= constants.min_fret;
        output(j,i+1) = median( data.fret(nonzero,j) );
    end

end %for each trace.

output(:,1) = output(:,1)/1000;  %convert to seconds


if nargout>0,
    varargout{1} = output;
    return;
end



%% Save the output to file.
if nFiles==1,
    [p,f] = fileparts(files{1});
    outputFilename = fullfile(p, [f '_' mfilename '.txt']);
else
    outputFilename = [mfilename '.txt'];
end

[f,p] = uiputfile('*.txt', [mfilename ': save output'], outputFilename);
outputFilename = fullfile(p,f);

if f~=0,
    % Write header line
    fid = fopen(outputFilename,'w');
    fprintf(fid,'Time (s)\t%s\n', strjoin(titles,'\t'));
    fclose(fid);

    dlmwrite(outputFilename,output,'-append','delimiter','\t');
end



%% Plot the result
ax = axes('Parent',figure);
plot( ax, output(:,1), output(:,2:end) );
xlabel(ax, 'Time (s)');
ylabel(ax, 'Average FRET value');
if nFiles>1,
    legend(ax, titles);
else
    title(ax, titles{1});
end




