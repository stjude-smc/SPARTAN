function [dwt,model] = vbFRET_dwt(filename)
% Converts a vbFRET idealization file (.mat) to a QuB dwell time file (DWT).

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Get filename from user
if nargin<1,
    [f,p] = uigetfile('*.mat','Select a vbFRET output file');
    if ~ischar(f)
        return;
    end
    filename = [p f];
end


% Load saved idealizations
vbIdl = load(filename);
paths = vbIdl.path;
FRET = vbIdl.FRET;
clear vbIdl;


% Load the data file that was used for idealization.
% This is only necessary for getting offsets and the sampling interval
% for the DWT file, which aren't saved by vbFRET.
[f p] = uigetfile('*.traces','Select original data file');
dataFile = [p f];

if ~ischar(f) || ~exist(dataFile,'file'),
    return;
end

data = loadTraces( dataFile );
traceLen = size(data.donor,2);
sampling = data.time(2)-data.time(1);
clear data;
% traceLen = 6000;
% sampling = 160;


% Iterate over traces, creating idealizations for each.
nTraces = numel( paths );
dwt = cell(nTraces,1);
model = cell(nTraces,1);

for i=1:nTraces,    
    path = paths{i};
    fretValues = unique( path );
    fretStd = repmat( 0.05, size(fretValues) ); %default value
    nFrames = numel(path);
    assert( nFrames<=traceLen, 'Idealized traces are longer than data!' );
    
    currentDwellLength = 1;
    stateAssignments = [];
    times = [];
    
    % Convert state assignments to state number (instead of FRET value).
    % Also, determine FRET standard deviations for states.
    for k=1:numel(fretValues),
        path( path==fretValues(k) ) = k;
        fretStd(k) = std( FRET{i}(path==k) );
    end
    fretStd = max( fretStd, repmat(0.01,size(fretStd)) ); %minimal std value
    
    % Iterate over frames within a trace, accumulating dwell time in the
    % current state until a transition is observed.
    for j=2:nFrames
        % A new dwell just started, add current dwell to DWT.
        if path(j) ~= path(j-1) || j==nFrames,
            times(end+1) = currentDwellLength;
            stateAssignments(end+1) = path(j-1);
            
            currentDwellLength = 1;
        else
            currentDwellLength = currentDwellLength+1;
        end
    end
    
    % Save results for this trace
    dwt{i} = [stateAssignments' times'];
    model{i} = [fretValues fretStd];
end

% Save as a dwell-time file
dwtFile = strrep( filename, '.mat', '.dwt' );
[f p] = uiputfile( '*.dwt', 'Save idealization file', dwtFile );

if ischar(f),
    dwtFile = [p f];
    offsets = ( 0:(nTraces-1) )*traceLen;
    saveDWT( dwtFile, dwt, offsets, model, sampling );
end



