function combineDatasets( filenames, outFilename, choice )
% combineDatasets  Combine several smFRET data files into one file.
%
%  combineDatasets( FILES, OUTPUT ) combines all of the smFRET data in the given
%  cell array of file names FILES, along any idealizations associated with these
%  files (same file name with a .qub.dwt extension) and saves the result as a
%  the new file OUTPUT. If not given, the user will be prompted for file names.
%
%  combineDatasets( ..., METHOD ) will combine traces files that are of
%  different lengths (number of frames). The options are 'truncate' to remove
%  the ends of longer traces or 'extend' to pad the ends of shorter traces so
%  that they all have the same length. If not given, the user will be prompted.

%   Copyright 2007-2015 Cornell University All Rights Reserved.



%% Get file names if not specified.
if nargin<1,
    filenames = getFiles;
end

nFiles = numel(filenames);
if nFiles<1, return; end



%% Get the data sizes to determine if the traces need to be resized.
[nTraces,traceLen] = sizeTraces(filenames);

if min(traceLen)~=max(traceLen),
    choices = {'Truncate','Extend','Cancel'};
    if nargin<3 || ~ismember(lower(choice),lower(choices)),
        choice = questdlg('Traces are different lengths. How do you want to modify the traces so they are all the same length?', ...
                  'Resize traces', choices{:}, 'Cancel');
    end
    
    if strcmpi(choice,'truncate'),
        newTraceLength = min(traceLen);
    elseif strcmpi(choice,'extend'),
        newTraceLength = max(traceLen);
    else
        return;  %user hit cancel
    end
else
    newTraceLength = traceLen(1);
end


%% Get output filename, if not specified.
if nargin<2 || isempty(outFilename),
    [f,p] = uiputfile('*.traces','Select output filename','combined.traces');
    if f==0, return; end
    outFilename = fullfile(p,f);
end


%% Load data
data = cell(nFiles,1);
idl  = cell(nFiles,1);
model= cell(nFiles,1);
sampling = zeros(nFiles,1);

h = waitbar(0,'Loading files...');

for i=1:nFiles,
    % Load traces from file.
    data{i} = loadTraces(filenames{i});
    
    % If it exists, load the associated dwell-time file (.dwt).
    [path,file] = fileparts( filenames{i} );
    dwt_fname = fullfile(path, [file '.qub.dwt']);
    
    if ~exist(dwt_fname,'file'),
        dwt_fname = fullfile(path, [file '.dwt']);
    end
    
    if exist(dwt_fname,'file'), 
        [dwt,sampling(i),offsets,model{i}] = loadDWT( dwt_fname );
        idl{i} = dwtToIdl(dwt, offsets, data{i}.nFrames, data{i}.nTraces);
        
        assert( sampling(i)==data{i}.sampling, 'Idealization time resolution mismatch' );
    else
        idl{i} = zeros(data{i}.nTraces, data{i}.nFrames);
    end
    
    waitbar(0.8*i/nFiles,h);
end


%% Resize traces so they are all the same length

waitbar(0.8,h,'Combining...');

if all(newTraceLength==traceLen),
    data = combine( data{:} );
    
else
    % Combine trace data, truncating or extending the data
    data = combine( data{:}, lower(choice) );
    
    % Truncate or pad idealization matrix to match the final size.
    for i=1:nFiles,
        delta = max(0, newTraceLength-traceLen(i) );
        padding = zeros(nTraces(i), delta);
        idl{i} = [ idl{i}(:,1:newTraceLength) padding ];
    end
end

% Combine idealizations.
idl = vertcat(idl{:});


%% Save merged traces and idealization data to file.

waitbar(0.9,h,'Saving...');

% Save merged trace data to file.
saveTraces( outFilename, data );


% Merge dwt files if present and consistent.    
nStates = cellfun( @numel, model ); %count number of states in each model.

if ~all( nStates(1)==nStates ),
    warning('Idealizations do not have the same number of states!');

elseif ~all( sampling(1)==sampling ),
    warning('Time resolution is not the same for all .dwt files!');

elseif ~all( idl(:)==0 )
    % Make an "average" model for the combined dwt file.
    modelAll = mean( cat(3,model{:}), 3 );  

    % Convert idealization to dwell-times for saving.
    [dwt,offsets] = idlToDwt(idl);

    % Save the dwt file.
    [p,f] = fileparts(outFilename);
    if isempty(p), p=pwd; end
    dwtFilename = fullfile(p, [f '.qub.dwt']);
    
    saveDWT( dwtFilename, dwt, offsets, modelAll, sampling(1) );
end


waitbar(1,h);
close(h);



