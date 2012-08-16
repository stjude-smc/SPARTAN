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
d = cell(0,1); a=d; f=d; ids=d; time=d;

for i=1:nFiles,

    % Load traces from file
    data = loadTraces(filenames{i});
    d{i} = data.donor;
    a{i} = data.acceptor;
    f{i} = data.fret;
    ids{i}  = data.ids;
    time{i} = data.time;
    
    assert( ~any(isnan(data.donor(:))) & ~any(isnan(data.acceptor(:))) & ~any(isnan(data.fret(:))) );
    
    nTraces = nTraces+size(data.donor,1);
    traceLen(i) = numel(data.time);
    assert( traceLen(i)>1 );
    
    waitbar(0.3*i/nFiles,h);
end
minTraceLen = min( traceLen );

% Resize traces so they are all the same length
for i=1:nFiles,
    if isempty(d{i}), continue; end
    d{i} = d{i}(:,1:minTraceLen);
    a{i} = a{i}(:,1:minTraceLen);
    f{i} = f{i}(:,1:minTraceLen);
    time{i} = time{i}(1:minTraceLen);
end



%%

d_out = [];
a_out = [];
f_out = [];
ids_out = {};

for i=1:nFiles,

    % Add data to combined dataset
    d_out = [d_out; d{i}];
    a_out = [a_out; a{i}];
    f_out = [f_out; f{i}];
    ids_out = [ids_out ; ids{i}];
    
    waitbar(0.3+0.2*i/nFiles,h);
end

assert( size(f_out,1)==nTraces );

clear data;
data.donor    = d_out;
data.acceptor = a_out;
data.fret     = f_out;
data.ids      = ids_out;
data.time     = time{1};

saveTraces( outFilename, 'traces', data );

waitbar(1,h);
close(h);
