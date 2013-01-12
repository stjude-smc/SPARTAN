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
    [f,p] = uiputfile('*.traces','Select output filename');
    if f==0, return; end
    outFilename = [p f];
end

nFiles = numel(filenames);
if nFiles<1, return; end


h = waitbar(0,'Combining datasets');


%% Load data
nTraces = 0;
traceLen = zeros(nFiles,1);
d = cell(0,1); a=d; f=d; time=[];
metadataAll = struct([]);

for i=1:nFiles,

    % Load traces from file
    data = loadTraces(filenames{i});
    d{i} = data.donor;
    a{i} = data.acceptor;
    f{i} = data.fret;
    
    % Find the longest time axis for merged traces.
    if length(data.time)>max(traceLen),
        time = data.time;
    end
    
    % Merge metadata fields into the final structure.
    if isfield(data,'traceMetadata'),
        if i==1,
            metadataAll = data.traceMetadata;
        else
            % Remove metadata fields that are not present in all datasets.
            % Otherwise, concatinating the two will give an error.
            data.traceMetadata = rmfield( data.traceMetadata, setdiff(fieldnames(data.traceMetadata),fieldnames(metadataAll)) );
            metadataAll   = rmfield( metadataAll, setdiff(fieldnames(metadataAll),fieldnames(data.traceMetadata)) );
            
            metadataAll = [metadataAll data.traceMetadata];
        end
    end
    
    assert( ~any(isnan(data.donor(:))) & ~any(isnan(data.acceptor(:))) & ~any(isnan(data.fret(:))) );
    
    nTraces = nTraces+size(data.donor,1);
    traceLen(i) = numel(data.time);
    assert( traceLen(i)>1 );
    
    waitbar(0.7*i/nFiles,h);
end




%% Resize traces so they are all the same length
if min( traceLen )~=max( traceLen ),
    
    choice = questdlg('Traces are different lengths. How do you want to modify the traces so they are all the same length?', ...
                      'Resize traces', 'Truncate','Extend','Cancel', 'Cancel');

    if strcmp(choice,'Cancel'),
        close(h);
        return;

    elseif strcmp(choice,'Truncate'),
        % Shorten traces to the same length
        newTraceLength = min( traceLen );

        for i=1:nFiles,
            if isempty(d{i}), continue; end
            d{i} = d{i}(:,1:newTraceLength);
            a{i} = a{i}(:,1:newTraceLength);
            f{i} = f{i}(:,1:newTraceLength);
        end
        
        time = time(1:newTraceLength);

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
        end
    end
    
end %if trace length mismatch

waitbar(0.8,h);


% Merge fluorescence and FRET data
data.time = time;
data.donor    = vertcat( d{:} );
data.acceptor = vertcat( a{:} );
data.fret     = vertcat( f{:} );
data.traceMetadata = metadataAll;

assert( size(data.fret,1)==nTraces );


% Save merged dataset to file.
[p,f,e] = fileparts(outFilename);
if ~isempty(strfind(e,'traces')),
    saveTraces( outFilename, 'traces', data );
elseif ~isempty(strfind(e,'txt')),
    saveTraces( outFilename, 'txt', data );
end

waitbar(1,h);
close(h);



