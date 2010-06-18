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

% TODO: 
%   - Rewrite the mem-map section.
%   * PROBLEM: the data cannot be written as it is loaded because of the way
%     it is (and must be) stored. One solution is to reserve space for the
%     data, memmap it, and then write it all in. To do so, we need to know
%     the size of the data all together - this can be estimated from the file
%     sizes. Create a larger file than anticipated so not expansion is
%     necessary. For .txt files, each number takes up ~10 bytes.
%   * PROBLEM: There is no easy way to store the IDs in the memory mapped data.
%     In order for a simple passthrough to work, the Data element must have
%     everything. As is, this is clearly impossible. What I can do is create
%     a more complex wrapper that will handle the text data differently.
%     
%   * Because this proceedure will be very slow, it should only be used in
%     extreme cases with LOTS of data. A second code path should be used in
%     other cases with all data loaded into memory.
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


% Load parameters
if nargin<2
    constants = cascadeConstants;
end

useMemmap = 0; %TEMP

wbh = waitbar(0,'Loading traces files...');

nFiles = numel(files);
nTraces = 0;
ids = {};
indexes = zeros( nFiles,2 ); %start/end trace indexes for each file.


if options.useMemmap,
    % First, estimate the space needed to store each file.
    nEl = filesize(files);
    totalSize = 1.2*sum(nEl); %with 20% safety margin.
    
    % Reserve space for the memory map, in groups of 1000 for speed.
    tempfile = [tempdir 'autotrace_memmap.dat'];
    fid = fopen(tempfile,'w');
    for i=1:ceil(totalSize/1000),
        fwrite(fid,zeros(1000,1),'double');
    end
    fclose(fid);
end


for i=1:numel(files),
    % Load traces data.
    [d,a,f,ids_t,time] = loadTraces( files{i} );
    
    % Verify traces files match in size.
    if exist('traceLen','var'),
        assert( traceLen==size(d,2), 'Trace size mismatch.' );
    else
        traceLen=size(d,2);
    end
    
    indexes(i,:) = [1 size(d,1)]+nTraces;
    
    
    % Save data
    if useMemmap,
        % Save data to file - will be loaded as a memory map later.
        fwrite( fid, d(1:nTraces,1:traceLen), 'double' );
        fwrite( fid, a(1:nTraces,1:traceLen), 'double' );
        fwrite( fid, f(1:nTraces,1:traceLen), 'double' );
    else
        % If this is the first run, allocate space for the remaining data, using
        % this file as a guess for the size.
        if i==1,
            d_all = zeros( size(d,1)*nFiles, size(d,2) );
            a_all = zeros( size(d,1)*nFiles, size(d,2) );
            f_all = zeros( size(d,1)*nFiles, size(d,2) );
        end
        
        d_all( nTraces+(1:size(d,1)), : ) = d;
        a_all( nTraces+(1:size(d,1)), : ) = a;
        f_all( nTraces+(1:size(d,1)), : ) = f;
    end
    
    nTraces = nTraces + size(d,1);
    ids = [ids ids_t];
    
    wbh = waitbar(0.9*i/nFiles,wbh);
end


if useMemmap,
    fclose(fid);
end




% Depending on performance settings, use a memory-mapped file instead of
% returning the data as one big variable in RAM.
if ~options.useMemmap
    data.d = d_all(1:nTraces,1:traceLen);
    data.a = a_all(1:nTraces,1:traceLen);
    data.f = f_all(1:nTraces,1:traceLen);
    data.time = time;
    data.ids = ids;
    
else
    % Load the data as a memmory-mapped variable.
    % The structure containing references to the memory mapped data is then
    % returned to the calling function. As far as I can tell, the data are not
    % hard copied into memory in this process...
    mmapdata = memmapfile( tempfile, 'Format', {...
                           'double',[nTraces traceLen],'d'; ...
                           'double',[nTraces traceLen],'a'; ...
                           'double',[nTraces traceLen],'f'; ...
                           'double',[1 traceLen],'time'; ...
                        }, 'Repeat',1,'Writable',false );

    data = mmapdata.Data;
end


% Clean up.
close(wbh);




end



function nEl = filesize( filenames )
% Returns a vector of the size (of each input file. The size is an estimate
% of the number of data elements in the data.


% Convert input to a cell array
if ~iscell( filesnames ),
    filenames = {filesnames};
end

nEl = zeros(1,numel(filenames));

for i=1:numel(filenames)
    d = dir( filenames{i} );
    assert( numel(d)==1, 'File does not exist' );
    
    % Estimate the number of elements
    [~,~,e] = fileparts( filenames{i} );
    
    if strcmp(e,'.txt')
        nEl(i) = d.bytes/10;
    elseif strcmp(e,'.traces')
        nEl(i) = d.bytes/2;
    else
        error('Unknown format - cannot estimate file size');
    end
end






end






