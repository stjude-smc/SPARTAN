function output = percentTime( filenames )
% PERCENTTIME  Stable state probabilities
%
%   Calculates the percentage time spent in each non-zero FRET state
%   using a QuB idealization data file (DWT).  States are listed in columns
%   from low to high FRET.  Each dataset is a row.
%   
% CAUTION: assumes first model state (1) is a zero-FRET state and removes it!


if nargin<1,
    
    filenames = cell(1,0);
    
    % Request filenames from user
    while 1,
        [file path]=uigetfile('*.dwt','Choose an idealization file:');
        if file==0, break; end  %user hit "cancel"

        filenames{end+1} = [path filesep file];
        disp( filenames{end} );
    end
end

nFiles = numel(filenames);

if nFiles == 0,
    disp('No files specified, exiting.');
    return;
end



output = zeros(0,0);

% Load dwell lifetime info from idealization data
for i=1:nFiles,
    
    totals = TotalTime( filenames{i} );
    
    % Sum up times over all molecules
    totals = totals(2:end);  %remove 0-FRET state
    
    assert( size(output,1)==0 || numel(totals) == size(output,2), ...
      'Model mismatch' );
    
    % Calculate and save percentages
    output(i,:) = 100*totals/sum(totals);
    
end %for each file





function totals = TotalTime( dwtfilename )

N = 0;
totals = zeros(1,4);

% Load dwell time file
fid = fopen(dwtfilename,'r');
assert( fid>0, 'ERROR: file does not exist!' );

while 1,
    % Load next segment in file  
    data = textscan(fid, 'Segment: %d %*[^\n]');
    segid = data{1};
    
    if numel(segid) == 0, break; end  %end of file
    N = N+1;
    
    data = textscan(fid, '%d%d');
    states = data{1}+1;
    times  = data{2};
    
    % Add times to total
    for i=1:numel(states),
        totals(states(i)) = totals(states(i)) + times(i);
    end
end

disp(N);




