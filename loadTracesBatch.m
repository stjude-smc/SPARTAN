function [data,ids,time] = loadTracesBatch( files )
% LOADBATCH   Loads a list of files or all traces in a directory.
%
%   [DATA,IDS,TIME] = loadBatch( FILES )
%   Loads fluorescence/FRET data from files in the cell array FILES. These
%   are combined together into a single binary array (DATA) that is stored
%   on disk and memory mapped. Memory mapping enables very fast access to
%   the data without having to store all the data in memory (the file
%   itself is in a temporary location).
%

%   [DATA,IDS,TIME] = loadBatch( DIR )
%   Load all .traces files in the directory DIR.
%

% Load parameters
constants = cascadeConstants;
options.useMemmap = constants.useMemmap;


nFiles = numel(files);
nTraces = 0;
ids = {};

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
    
    % Save data
    d_all( nTraces+(1:size(d,1)), : ) = d;
    a_all( nTraces+(1:size(d,1)), : ) = a;
    f_all( nTraces+(1:size(d,1)), : ) = f;
    nTraces = nTraces + size(d,1);
    
    ids = [ids ids_t];
end


% Depending on performance settings, use a memory-mapped file instead of
% returning the data as one big variable in RAM.
if ~options.useMemmap
    % This will use up more RAM, but is faster because the data do
    % not have to be saved to file.
    data.d = d_all(1:nTraces,1:traceLen);
    data.a = a_all(1:nTraces,1:traceLen);
    data.f = f_all(1:nTraces,1:traceLen);

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




