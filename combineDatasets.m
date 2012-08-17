function combineDatasets( filenames, outFilename )
% combineDatasets  Combine several smFRET data files into one file.
%
%  combineDatasets( FILES, OUTPUT )
%  Each smFRET data file (see loadTraces.m) specified in the the cell array
%  FILES is loaded and combined into a single, large dataset. The resulting
%  dataset is saved to the OUTPUT filename. If files are of differing
%  lengths, they are truncated to the minimal size. If input arguments are
%  not specified, the user will be prompted for them.


%%
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
d = cell(0,1); a=d; f=d; time=d;

for i=1:nFiles,

    % Load traces from file
    data = loadTraces(filenames{i});
    d{i} = data.donor;
    a{i} = data.acceptor;
    f{i} = data.fret;
    time = data.time;
    
    if i==1,
        metadataAll = data.traceMetadata;
    else
        metadataAll = [metadataAll data.traceMetadata];
    end
    
    assert( ~any(isnan(data.donor(:))) & ~any(isnan(data.acceptor(:))) & ~any(isnan(data.fret(:))) );
    
    nTraces = nTraces+size(data.donor,1);
    traceLen(i) = numel(data.time);
    assert( traceLen(i)>1 );
    
    waitbar(0.7*i/nFiles,h);
end
minTraceLen = min( traceLen );


% Resize traces so they are all the same length
for i=1:nFiles,
    if isempty(d{i}), continue; end
    d{i} = d{i}(:,1:minTraceLen);
    a{i} = a{i}(:,1:minTraceLen);
    f{i} = f{i}(:,1:minTraceLen);
end
waitbar(0.8,h);


% Merge fluorescence and FRET data
data.time = time(1:minTraceLen);
data.donor    = vertcat( d{:} );
data.acceptor = vertcat( a{:} );
data.fret     = vertcat( f{:} );
data.traceMetadata = metadataAll;

assert( size(data.fret,1)==nTraces );


% Save merged dataset to file.
saveTraces( outFilename, 'traces', data );

waitbar(1,h);
close(h);



