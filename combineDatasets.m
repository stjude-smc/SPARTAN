function output = combineDatasets( filenames, outFilename )
% combineDatasets  Combine several smFRET data files into one file.
%
%  combineDatasets( FILES, OUTPUT )
%  Each smFRET data file (see loadTraces.m) specified in the the cell array
%  FILES is loaded and combined into a single, large dataset. The resulting
%  dataset is saved to the OUTPUT filename. If files are of differing
%  lengths, they are truncated to the minimal size. If input arguments are
%  not specified, the user will be prompted for them.
%


%% Get file names if not specified.
if nargin<1,
    filenames = getFiles;
    if isempty(filenames), return; end;
end

if nargin<2,
    [f,p] = uiputfile('*.traces','Select output filename','combined.traces');
    if f==0, return; end
    outFilename = [p f];
end

nFiles = numel(filenames);
if nFiles<1, return; end


h = waitbar(0,'Combining datasets');


%% Load data
nTraces  = zeros(1,nFiles);
traceLen = zeros(1,nFiles);

data = cell(nFiles,1);
dwt=cell(nFiles,1);
sampling=zeros(nFiles,1);
offsets=cell(nFiles,1);
model=cell(nFiles,1);

for i=1:nFiles,

    % Load traces from file. If any of the files have a slightly different
    % structure (fields missing or in a different order), this will fail!
    data{i} = loadTraces(filenames{i});
    [nTraces(i),traceLen(i)] = size(data{i}.donor);
    
    % Look for an idealization. These will be combined if every dataset has one.
    [path,file] = fileparts( filenames{i} );
    dwt_fname = [path filesep file '.qub.dwt'];
    if exist(dwt_fname,'file'), 
        [dwt{i},sampling(i),offsets{i},model{i}] = loadDWT( dwt_fname );
    end
    
    waitbar(0.7*i/nFiles,h);
end


%% Resize traces so they are all the same length
newTraceLength = traceLen(1);

if min(traceLen) ~= max(traceLen),
    
    choice = questdlg('Traces are different lengths. How do you want to modify the traces so they are all the same length?', ...
                  'Resize traces', 'Truncate','Extend','Cancel', 'Cancel');

    % ---- User hit cancel, do nothing.
    if strcmp(choice,'Cancel'),
        close(h);
        return;
        
    % ---- Truncate traces to the same length
    elseif strcmp(choice,'Truncate'),
        newTraceLength = min(traceLen);
        
        output = combine( data{:} );

        for i=1:nFiles,
            % Convert dwell-times to an idealization, which are easy to
            % truncate, and convert back for saving later.
            if ~isempty(dwt{i}),
                idl = dwtToIdl( dwt{i}, traceLen(i), offsets{i} );
                idl = idl(:,1:newTraceLength);
                [dwt{i},offsets{i}] = idlToDwt(idl);
            end
        end %end for each file.

        
    % ---- Extend traces to the same length
    else
        % Extend the traces so they are the same length by duplicating
        % the values in the last frame to pad the end. This can create 
        % problems later on, so use with caution!
        newTraceLength = max( traceLen );
        
        output = combine( data{:},'extend' );

        for i=1:nFiles,
            % To "extend" the dwell-times, we just have to change the
            % offsets. But this is very easy to do by converting it to an
            % idealization and adding zeros (a marker for regions that are
            % not idealized) to the ends.
            if ~isempty(dwt{i}),
                delta = newTraceLength - traceLen(i);
                
                idl = dwtToIdl( dwt{i}, traceLen(i), offsets{i} );
                idl = [ idl  zeros( size(idl,1), delta )  ];
                [dwt{i},offsets{i}] = idlToDwt(idl);
            end
        end
    end
    
else
    output = combine( data{:} );
end %if trace length mismatch


clear data;
waitbar(0.9,h);


% Save merged dataset to file.
[p,f,e] = fileparts(outFilename);
if ~isempty(strfind(e,'traces')),
    saveTraces( outFilename, 'traces', output );
elseif ~isempty(strfind(e,'txt')),
    saveTraces( outFilename, 'txt', output );
end


% Merge dwt files if present and consistent.    
n = cellfun( @numel, model ); %count number of states in each model.
if ~all( n(1)==n ),
    warning('.dwt files found for all files, but models have different numbers of states!');

elseif ~all( sampling(1)==sampling ),
    warning('.dwt files found for all files, but they are not the same time resolution and cannot be combined.');

elseif ~all( cellfun(@isempty,dwt) )
    disp('Combining dwt files. I hope you used the same models for these!');

    offsetsAll = [];
    dwtAll     = {};
    modelAll   = zeros( size(model{1}) );

    % Merge all the idealizations, adjusting the offsets.
    fileOffsets = cumsum( [0 nTraces.*newTraceLength] );
    for i=1:numel(offsets),
        offsetsAll = [offsetsAll offsets{i}+fileOffsets(i) ];
        dwtAll     = [dwtAll dwt{i}];
        modelAll   = modelAll + model{i};
    end
    modelAll = modelAll./numel(offsets); %this gives us an "average" model.

    % Save the dwt file.
    if isempty(p), p=pwd; end
    dwtFilename = [p filesep f '.qub.dwt'];
    saveDWT( dwtFilename, dwtAll, offsetsAll, modelAll, sampling(1) );
end


waitbar(1,h);
close(h);



