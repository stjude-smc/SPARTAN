function output = avgFretTime( files )
% Averages the FRET values of all traces at each time point to generate a
% timecourse of changes in FRET over the population. Occupancy in zero-FRET
% states are not added into the sum, so photobleaching doesn't interfere.
% This is useful for stopped-flow drug/ligand delivery to measure the rate
% of change in the distribution. Also see stateTimecourse.m.
%


% Settings
truncateLen = 100;  %frames to calculate over
constants.min_fret = 0.2;  % minimum fret value, below which we assume there is no FRET.


% Get list of files if not specified.
if nargin<1 || isempty(files),
    files = getFiles;
end

if ~iscell(files), files = {files}; end

nFiles = numel(files);


% 
output = zeros(truncateLen,1+nFiles);

for i=1:numel(files),
    % Extra FRET traces from the current file.
    data = loadTraces( files{i} );
    f = data.fret(:,1:truncateLen);
    
    if i==1,
        output(:,1) = data.time(1:truncateLen);
    else
        if ~all(data.time(1:truncateLen)==output(:,1)),
            warning( 'Time axes do not match!' );
        end
    end
    
    % For each trace, average the FRET values at each frame to create an
    % average FRET trace. Exclude regions where the dyes are dark.
    for j=1:truncateLen,
        nonzero = f(:,j) >= constants.min_fret;
        output(j,i+1) = mean( f(nonzero,j) );
    end

end %for each trace.

output(:,1)=output(:,1)/1000; %convert to seconds



% Save the output to file.
if nFiles==1,
    [p,f] = fileparts(files{1});
    outputFilename = [p filesep f '_avgFretTime.txt'];
else
    [f,p] = uiputfile('*.txt','Save state timecourse data','avgFretTime.txt');
    outputFilename = [p f];
end

fid = fopen(outputFilename,'w');

for n=1:nFiles,
    if nFiles == 1,
        f = 'Mean FRET';
    else
        [~,f] = fileparts(files{n});
    end
    
    if output(1,1)==1,
        fprintf(fid,'Time (Frames)\t%s',f);
    else
        fprintf(fid,'Time (s)\t%s',f);
    end
end
fprintf(fid,'\n');
fclose(fid);

dlmwrite(outputFilename,output,'-append','delimiter','\t');



