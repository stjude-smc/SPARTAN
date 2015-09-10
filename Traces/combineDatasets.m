function output = combineDatasets( filenames, outFilename, choice )
% combineDatasets  Combine several smFRET data files into one file.
%
%  combineDatasets( FILES, OUTPUT )
%  Each smFRET data file (see loadTraces.m) specified in the the cell array
%  FILES is loaded and combined into a single, large dataset. The resulting
%  dataset is saved to the OUTPUT filename. If files are of differing
%  lengths, they are truncated to the minimal size. If input arguments are
%  not specified, the user will be prompted for them.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


output = [];


%% Get file names if not specified.
if nargin<1,
    filenames = getFiles;
    if isempty(filenames), return; end;
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
        return;
    end
else
    newTraceLength = traceLen(1);
end


% Get output filename, if not specified.
if nargin<2 || isempty(outFilename),
    [f,p] = uiputfile('*.traces','Select output filename','combined.traces');
    if f==0, return; end
    outFilename = [p f];
end


%% Load data
data = cell(nFiles,1);
dwt=cell(nFiles,1);
sampling=zeros(nFiles,1);
offsets=cell(nFiles,1);
model=cell(nFiles,1);

h = waitbar(0,'Combining datasets');

for i=1:nFiles,
    % Load traces from file.
    data{i} = loadTraces(filenames{i});
    
    % Look for an idealization. These will be combined if every dataset has one.
    [path,file] = fileparts( filenames{i} );
    dwt_fname = fullfile(path, [file '.qub.dwt']);
    if exist(dwt_fname,'file'), 
        [dwt{i},sampling(i),offsets{i},model{i}] = loadDWT( dwt_fname );
    end
    
    waitbar(0.7*i/nFiles,h);
end


%% Resize traces so they are all the same length

if all(newTraceLength==traceLen),
    % Combine trace data
    output = combine( data{:} );
    
else
    % Combine trace data, truncating or extending the data
    output = combine( data{:}, lower(choice) );
    
    % Truncate or extend idealizations.
    % FIXME: might be useful as a separate function
    for i=1:nFiles,
        if isempty(dwt{i}), continue; end
        
        % Convert dwell-times to an idealization, which are easy to
        % truncate, and convert back for saving later.
        idl = dwtToIdl( dwt{i}, traceLen(i), offsets{i} );

        % This will extend (with zeros), if applicable.
        delta = max(0, newTraceLength-traceLen(i) );
        idl = [ idl  zeros( size(idl,1), delta )  ];  %#ok

        % This will truncate, if applicable.
        idl = idl(:,1:newTraceLength);

        [dwt{i},offsets{i}] = idlToDwt(idl);
    end
end

waitbar(0.9,h);


% Save merged trace data to file.
saveTraces( outFilename, output );


% Merge dwt files if present and consistent.    
n = cellfun( @numel, model ); %count number of states in each model.
if ~all( n(1)==n ),
    warning('.dwt files found for all files, but models have different numbers of states!');

elseif ~all( sampling(1)==sampling ),
    warning('.dwt files found for all files, but they are not the same time resolution and cannot be combined.');

elseif ~all( cellfun(@isempty,dwt) )
    disp('Combining dwt files. I hope you used the same models for these!');

    % Make an "average" model for the combined dwt file.
    modelAll = mean( cat(3,model{:}), 3 );  

    % Update offsets for combined file.
    fileOffsets = cumsum( [0; nTraces*newTraceLength] );
    for i=1:nFiles,
        offsets{i} = offsets{i} + fileOffsets(i);
    end

    % Save the dwt file.
    [p,f] = fileparts(outFilename);
    if isempty(p), p=pwd; end
    dwtFilename = fullfile(p, [f '.qub.dwt']);
    saveDWT( dwtFilename, [dwt{:}], [offsets{:}], modelAll, sampling(1) );
end


waitbar(1,h);
close(h);



