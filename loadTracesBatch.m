function [data,indexes] = loadTracesBatch( files, constants )
% LOADBATCH   Loads a list of files or all traces in a directory.
%
%   [DATA,INDEXES] = loadBatch( FILES )
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
    if isempty(files), return; end
end

% Load parameters
if nargin<2
    constants = cascadeConstants;
end

options.useMemmap = constants.useMemmap;

wbh = waitbar(0,'Loading traces files...');

nFiles = numel(files);
nTraces = 0;
ids = {};
indexes = zeros( nFiles,2 ); %start/end trace indexes for each file.

for i=1:numel(files),
    % Load traces data.
    [d,a,f,ids_t,time] = loadTraces( files{i} );
    
    % Verify traces files match in size.
    if exist('traceLen','var'),
        assert( traceLen==size(d,2), 'Trace size mismatch.' );
    else
        traceLen=size(d,2);
    end
    
    % If this is the first run, allocate space for the remaining data, using
    % this file as a guess for the size.
    if i==1,
        d_all = zeros( size(d,1)*nFiles, size(d,2) );
        a_all = zeros( size(d,1)*nFiles, size(d,2) );
        f_all = zeros( size(d,1)*nFiles, size(d,2) );
    end
    
    indexes(i,:) = [1 size(d,1)]+nTraces;
    
    
    % Save data
    d_all( nTraces+(1:size(d,1)), : ) = d;
    a_all( nTraces+(1:size(d,1)), : ) = a;
    f_all( nTraces+(1:size(d,1)), : ) = f;
    nTraces = nTraces + size(d,1);
    
    ids = [ids ids_t];
    
    wbh = waitbar(0.9*i/nFiles,wbh);
end


% Depending on performance settings, use a memory-mapped file instead of
% returning the data as one big variable in RAM.
if ~options.useMemmap
    data.d = d_all(1:nTraces,1:traceLen);
    data.a = a_all(1:nTraces,1:traceLen);
    data.f = f_all(1:nTraces,1:traceLen);
    data.ids = ids;
    data.time = time;
    
else
    tempfile = [tempdir 'autotrace_memmap.dat'];
    
    % Write binary data to file:
    fid = fopen(tempfile,'w');
    fwrite( fid, d_all(1:nTraces,1:traceLen), 'double' );
    fwrite( fid, a_all(1:nTraces,1:traceLen), 'double' );
    fwrite( fid, f_all(1:nTraces,1:traceLen), 'double' );
    fclose(fid);

    % Load the data as a memmory-mapped variable.
    % The structure containing references to the memory mapped data is then
    % returned to the calling function. As far as I can tell, the data are not
    % hard copied into memory in this process...
    mmapdata = memmapfile( tempfile, 'Format', {...
                           'double',[nTraces traceLen],'d'; ...
                           'double',[nTraces traceLen],'a'; ...
                           'double',[nTraces traceLen],'f'; ...
                        }, 'Repeat',1,'Writable',false );

    data = mmapdata.Data;
end


% Clean up.
close(wbh);


