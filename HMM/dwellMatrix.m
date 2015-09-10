function dwells = dwellMatrix(dwtfilename)
% DWELL_MATRIX  Loads a QuB idealization dwell time file
% 
%    Asks user for DWT filename to load and outputs a NxM matrix of N
%    molecules and M states containing the total amount of time the
%    molecule spent in that state (in ms).
%    Outfile filename will be xxxxx.dwt.txt

%   Copyright 2007-2015 Cornell University All Rights Reserved.


if ~exist( 'dwtfilename', 'var' ),
    [datafile,datapath] = uigetfile({'*.dwt'},'Choose a DWT file:');
    if datafile==0, return; end  %user hit "cancel"

    dwtfilename = fullfile(datapath,datafile);
end


dwells = zeros(1,1);  %will expand automatically with zeros

% Load dwell time file
fid = fopen(dwtfilename,'r');
assert( fid>0, 'ERROR: file does not exist!' );

while 1,
    % Load next segment in file  
    data = textscan(fid, 'Segment: %d %*[^\n]');
    segid = data{1};
    
    if numel(segid) == 0, break; end  %end of file
    
    data = textscan(fid, '%d%d');
    states = data{1}+1;
    times  = data{2};
    
    % Expand matrix for next row
    dwells(segid,max(states)) = 0;
    
    % Save dwells in segment
    for i=1:numel(states),
        dwells( segid, states(i) ) = dwells( segid, states(i) ) + double(times(i)');
    end
end
fclose(fid);


% Write matrix to file
%out_filename = strrep( dwtfilename, '.dwt', '.dwt.txt' );

%dlmwrite(out_filename, dwells, ' ');


end  % function DWELL_MATRIX
