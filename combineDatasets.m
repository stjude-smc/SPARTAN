function combineDatasets( filenames, outFilename )

%%
if nargin<1,
    filenames = getFiles('*.txt');
end

if nargin<2,
    [f,p] = uiputfile('*.txt','Select output filename');
    if f==0, return; end
    outFilename = [p f];
end

nFiles = numel(filenames);

if nFiles<1, return; end


%% Load data
nTraces = 0;
traceLen = zeros(nFiles,1);
d = cell(0,1); a=d; f=d; ids=d; time=d;

for i=1:nFiles,

    % Load traces from file
    [d_in,a_in,f_in,ids_in,time_in] = loadTraces(filenames{i});
    d{i} = d_in;
    a{i} = a_in;
    f{i} = f_in;
    ids{i} = ids_in;
    time{i} = time_in;
    
    nTraces = nTraces+size(d{i},1);
    traceLen(i) = size(d{i},2);
    
    if traceLen(i)<1,
        error( ['File is empty: ' filenames{i}] );
    end
end
minTraceLen = min( traceLen );

% Resize traces so they are all the same length
for i=1:nFiles,
    d{i} = d{i}(:,1:minTraceLen);
    a{i} = a{i}(:,1:minTraceLen);
    f{i} = f{i}(:,1:minTraceLen);
    time{i} = time{i}(:,1:minTraceLen);
end



%%

d_out = [];
a_out = [];
f_out = [];
ids_out = {};

h = waitbar(0,'Combining datasets');

for i=1:nFiles,

    % Add data to combined dataset
    d_out = [d_out; d{i}];
    a_out = [a_out; a{i}];
    f_out = [f_out; f{i}];
    ids_out = [ids_out ids{i}];
    
    waitbar(0.9*i/nFiles,h);
end

assert( size(f_out,1)==nTraces );
saveTraces( outFilename, 'txt', d_out,a_out,f_out,ids_out,time{1} );

forQuB2( {outFilename} );

waitbar(1,h);
close(h);
