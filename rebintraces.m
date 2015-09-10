function rebintraces(factor, files)
%REBINTRACES  Creates traces with (simulated) reduced time resolution
% 
%   REBINTRACES(FILES)
%   Shrinks the temporal resolution of trace data by FACTOR (0.1 means half
%   the original resolution)
%   If no FILES are specified, user will be asked to select them.

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% TODO: add varargin options for stuff


% INITIALIZE & PROCESS FUNCTION PARAMETERS
assert( factor>=1, 'Cannot expand trace resolution' );


%% INITIALIZE & PROCESS FUNCTION PARAMETERS

% If no files specified, prompt user for them.
if ~exist('files','var'),
    disp('Select traces files, hit cancel when finished');
    
    filter = {'*.traces;*.rawtraces','Binary Traces Files (*.traces,*.rawtraces)'; ...
          '*.txt','Text Files (*.txt)'; ...
          '*.*','All Files (*.*)'};
    
    files = getFiles(filter);
end

nFiles = numel(files);

if nFiles == 0,
    disp('No files specified, exiting.');
    return;
end


%% 
for i=1:nFiles,
    
    % Load traces file
    data = LoadTraces( files{i} );
    
    % Shrink it
    newSize = ceil(size(data.donor,2)/factor);
    nTraces = size(data.donor,1);
    disp( sprintf('Resizing to %d frames',newSize) );
    
    data2.donor = zeros(nTraces,newSize);
    data2.acceptor = zeros(nTraces,newSize);
    
    for j=1:nTraces,
       data2.donor(j,:)    = TimeScale( data.donor(j,:),    factor );
       data2.acceptor(j,:) = TimeScale( data.acceptor(j,:), factor );  
    end
    data2.times = times(1:newSize); %*factor;
    
    %Should actually have correctTraces here??
    data2.fret = data2.acceptor./(data2.donor+data2.acceptor);
    
    
    % Write file back to disk
    [p f] = fileparts( files{i} );
    saveTraces( fullfile(p, [f '_rebin.traces']), data2 );
end

% end function


function output = TimeScale( trace, factor )

assert( factor >= 1 );

nFrames2 = ceil(numel(trace)/factor);
output = zeros( 1,nFrames2 );

j = 1;
itr2 = 1;
while floor(j)+1 <= numel(trace),
    
    itr1 = floor(j);
    
    r1 = 1 - (j-floor(j));
    r2 = factor-r1;
    
    output(itr2) = (r1*trace(itr1) + r2*trace(itr1+1))/factor;
    
    j = j+factor;
    itr2 = itr2+1;
    
end










