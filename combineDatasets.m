function combineDatasets( filenames, outFilename )
% combineDatasets  Combine several smFRET data files into one file.
%
%  combineDatasets( FILES, OUTPUT )
%  Each smFRET data file (see loadTraces.m) specified in the the cell array
%  FILES is loaded and combined into a single, large dataset. The resulting
%  dataset is saved to the OUTPUT filename. If files are of differing
%  lengths, they are truncated to the minimal size. If input arguments are
%  not specified, the user will be prompted for them.



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
nTracesTotal = 0;
nTraces  = zeros(1,nFiles);
traceLen = zeros(1,nFiles);
d = cell(0,1); a=d; f=d; time=[];
metadataAll = struct([]);
dwt=cell(nFiles,1); sampling=zeros(nFiles,1);
offsets=cell(nFiles,1); model=cell(nFiles,1);

for i=1:nFiles,

    % Load traces from file
    data = loadTraces(filenames{i});
    d{i} = data.donor;
    a{i} = data.acceptor;
    f{i} = data.fret;
    
    % Look for an idealization. These will be combined if every dataset has one.
    [path,file] = fileparts( filenames{i} );
    dwt_fname = [path filesep file '.qub.dwt'];
    if exist(dwt_fname,'file'), 
        [dwt{i},sampling(i),offsets{i},model{i}] = loadDWT( dwt_fname );
    end
    
    % Find the longest time axis for merged traces.
    if length(data.time)>max(traceLen),
        time = data.time;
    end
    
    % Merge metadata fields into the final structure.
    assert( isfield(data,'traceMetadata'), 'File doesn''t have metadata. This should never happen!' );
    if i==1,
        metadataAll = data.traceMetadata;
    else
        % Remove metadata fields that are not present in all datasets.
        % Otherwise, concatinating the two will give an error.
        data.traceMetadata = rmfield( data.traceMetadata, setdiff(fieldnames(data.traceMetadata),fieldnames(metadataAll)) );
        metadataAll        = rmfield( metadataAll, setdiff(fieldnames(metadataAll),fieldnames(data.traceMetadata)) );

        metadataAll = [metadataAll data.traceMetadata];
    end
    
    assert( ~any(isnan(data.donor(:))) & ~any(isnan(data.acceptor(:))) & ~any(isnan(data.fret(:))) );
    
    nTraces(i) = size(data.donor,1);
    nTracesTotal = nTracesTotal+size(data.donor,1);
    traceLen(i) = numel(data.time);
    assert( traceLen(i)>1 );
    
    waitbar(0.7*i/nFiles,h);
end




%% Resize traces so they are all the same length
newTraceLength = traceLen(1);

if min( traceLen )~=max( traceLen ),
    
    choice = questdlg('Traces are different lengths. How do you want to modify the traces so they are all the same length?', ...
                      'Resize traces', 'Truncate','Extend','Cancel', 'Cancel');

    % ---- User hit cancel, do nothing.
    if strcmp(choice,'Cancel'),
        close(h);
        return;
        
    % ---- Truncate traces to the same length
    elseif strcmp(choice,'Truncate'),
        
        newTraceLength = min( traceLen );

        for i=1:nFiles,
            if isempty(d{i}), continue; end
            d{i} = d{i}(:,1:newTraceLength);
            a{i} = a{i}(:,1:newTraceLength);
            f{i} = f{i}(:,1:newTraceLength);
            
            % Convert dwell-times to an idealization, which are easy to
            % truncate, and convert back for saving later.
            if ~isempty(dwt{i}),
                idl = dwtToIdl( dwt{i}, traceLen(i), offsets{i} );
                idl = idl(:,1:newTraceLength);
                [dwt{i},offsets{i}] = idlToDwt(idl);
            end
        end %end for each file.
        
        time = time(1:newTraceLength);

        
    % ---- Extend traces to the same length
    else
        % Extend the traces so they are the same length by duplicating
        % the values in the last frame to pad the end. This can create 
        % problems later on, so use with caution!
        newTraceLength = max( traceLen );

        for i=1:nFiles,
            if isempty(d{i}), continue; end

            delta = newTraceLength-size(d{i},2);
            d{i} = [d{i} repmat( d{i}(:,end), 1, delta )];
            a{i} = [a{i} repmat( a{i}(:,end), 1, delta )];
            f{i} = [f{i} zeros( size(d{i},1), delta )];
        
            % To "extend" the dwell-times, we just have to change the
            % offsets. But this is very easy to do by converting it to an
            % idealization and adding zeros (a marker for regions that are not
            % idealized) to the ends.
            if ~isempty(dwt{i}),
                idl = dwtToIdl( dwt{i}, traceLen(i), offsets{i} );
                idl = [ idl  zeros( size(idl,1), delta )  ];
                [dwt{i},offsets{i}] = idlToDwt(idl);
            end
        end
    end
    
end %if trace length mismatch

waitbar(0.75,h);


% Merge fluorescence and FRET data
data.time = time;
data.donor    = vertcat( d{:} );
data.acceptor = vertcat( a{:} );
data.fret     = vertcat( f{:} );
data.traceMetadata = metadataAll;

assert( size(data.fret,1)==nTracesTotal );
waitbar(0.9,h);


% Save merged dataset to file.
[p,f,e] = fileparts(outFilename);
if ~isempty(strfind(e,'traces')),
    saveTraces( outFilename, 'traces', data );
elseif ~isempty(strfind(e,'txt')),
    saveTraces( outFilename, 'txt', data );
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



