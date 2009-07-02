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

%%

d_out = [];
a_out = [];
f_out = [];
ids_out = {};

h = waitbar(0,'Combining datasets');

nTraces = 0;

for i=1:nFiles,

    % Load traces from file
    [d,a,f,ids,time] = loadTraces(filenames{i});
    nTraces = nTraces+size(d,1);

    % Add data to combined dataset
    d_out = [d_out; d];
    a_out = [a_out; a];
    f_out = [f_out; f];
    ids_out = [ids_out ids];
    
    waitbar(0.9*i/nFiles,h);
end

assert( size(f_out,1)==nTraces );
saveTraces( outFilename, 'txt', d_out,a_out,f_out,ids_out,time );

forQuB2( {outFilename} );

waitbar(1,h);
close(h);
