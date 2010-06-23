function [data,ids,indexes] = loadTracesBatch( files )
% LOADBATCH   Loads a list of files or all traces in a directory.
%
%   [DATA,IDS,INDEXES] = loadBatch( FILES )
%   Loads fluorescence/FRET data from files in the cell array FILES. These
%   are combined together into a structure (DATA) containing binary arrays
%   for donor, acceptor, and fret (d,a,f) data; the time axis (.time) and
%   a cell array of ids (.ids). INDEXES is an Nx2 array of the start and end
%   trace IDs (in the large arrays) for each file. This allows one to easily
%   seperate out the traces contributed by each file in the list.
%

%   [DATA,IDS,TIME] = loadBatch( DIR )
%   Load all .traces files in the directory DIR.
%


if nargin<1,
    filter = {'*.txt;*.traces','All Traces Files (*.txt,*.traces)'; ...
              '*.txt','Text Files (*.txt)'; ...
              '*.traces','Binary Traces Files (*.traces)'; ...
              '*.*','All Files (*.*)'};
    files = getFiles( filter );
    if isempty(files),
        data = [];
        return;
    end
end


wbh = waitbar(0,'Loading traces files...');
nFiles = numel(files);


% 1. Load up the trace headers so we know how much space needs to be allocated.
nTraces  = zeros( nFiles,1 );
traceLen = 0;

for i=1:nFiles,
    traceLen_last = traceLen;
    [nTraces(i),traceLen] = traceFileSize( files{i} );

    % Check that trace lenghts all match.
    if traceLen_last~=0 && traceLen~=traceLen_last,
        error('loadTraces: all traces must be of the same length!');
    end
end

nTracesTotal = sum(nTraces);


% Use memory mapping only if loading a very large amount of data
% useMemmap = (nTracesTotal > 5000);
useMemmap = 1;


% 2. Determine total data size and allocate space for all the trace data.
if useMemmap,
    disp('Using memory mapping to load large dataset.');
    
    % Reserve space for the memory map, in groups of 1000 for speed.
    tempfile = [tempdir 'autotrace_memmap.dat'];
    fid = fopen(tempfile,'w');
    
    totalSize = (1+3*nTracesTotal)*traceLen;
    for i=1:ceil(totalSize/1000),
        fwrite(fid,zeros(1000,1),'double');
    end
    fclose(fid);
    
    % Load the memory map:
    mmap = memmapfile( tempfile, 'Format', {...
                           'double',[nTracesTotal traceLen],'d'; ...
                           'double',[nTracesTotal traceLen],'a'; ...
                           'double',[nTracesTotal traceLen],'f'; ...
                           'double',[1 traceLen],'time'; ...
                        }, 'Repeat',1,'Writable',true );
    data = mmapPassthrough( mmap );
    
    
else
    data.d = zeros( nTracesTotal, traceLen );
    data.a = zeros( nTracesTotal, traceLen );
    data.f = zeros( nTracesTotal, traceLen );
end

ids = cell( nTracesTotal,1 );



% 3. Load up the data into memory.
nTraces = 0;
indexes = zeros( nFiles,2 ); %start/end trace indexes for each file.

for i=1:numel(files),
    % Load traces data.
    [d,a,f,ids_t,time] = loadTraces( files{i} );
        
    indexes(i,:) = [1 size(d,1)]+nTraces;
    
    % Verify time axes match
    %...
    
    % Save data
    idx = nTraces+(1:size(d,1));
    data.d( idx, : ) = d;
    data.a( idx, : ) = a;
    data.f( idx, : ) = f;
    ids( idx ) = ids_t;
    
    nTraces = nTraces + size(d,1);
    
    wbh = waitbar(0.9*i/nFiles,wbh);
end

data.time = time;



% Clean up.
close(wbh);

end



function [nTraces,traceLen] = traceFileSize( filename )
% Returns,,,

fid = fopen( filename, 'r' );
traceLen=fread(fid,1,'int32');
nTraces=fread(fid,1,'int16');
fclose(fid);

end









