function saveTraces( filename, format, varargin )
% SAVETRACES    Saves trace data to file
%
%    SAVETRACES( FNAME, 'txt',    D,A,F,IDs, time )
%    SAVETRACES( FNAME, 'traces', D,A,F,IDs, time, metadata )
%    SAVETRACES( FNAME, 'qub',    F )
%
% D:    MxN Donor fluorescence intensity matrix
% A:    MxN Acceptor ...
% F:    MxN FRET matrix
% IDs:  Mx1 cell array (strings) of trace identifiers
%      (M=traces,N=datapoints)
% Time:
% 

Nargs = numel(varargin);

switch format
    case 'txt'
        if Nargs<3, error('SAVETRACES: required parameters missing'); end
        saveTracesTxt( filename, varargin{:} );
        
    case 'traces'
        if Nargs<2, error('SAVETRACES: required parameters missing'); end
        saveTracesBinary( filename, varargin{:} );
        
    case 'qub'
        if Nargs~=1, error('SAVETRACES: required parameters missing'); end
        saveTracesQUB( filename, varargin{:} );
end

% END FUNCTION SAVETRACES



function saveTracesTxt( filename, donor,acceptor,fret, ids, time )
% FORMAT:
%   1 2 3 4 ... N
%   ID1 donor data
%   ID1 acceptor data
%   ID1 FRET data
%   ID2 donor data
%   ...
% 

[Ntraces,tlen] = size(donor);

if ~exist('time','var'),
    time = 1:tlen; %in frames
end

% Create IDs if not specified
if ~exist('ids','var'),
    [p,name] = fileparts(filename);
    
    ids = cell(Ntraces,1);
    for j=1:Ntraces;
        ids{j} = sprintf('%s_%d', name, j);
    end
end

% Remove special characters from IDs
ids = strrep( ids, '-', '_' );      %- is used as ID seperator, must be avoided
ids = strrep( ids, ' ', '_' );      %auto.txt format doesn't allow spaces

% Verify input arguments
if any( size(donor)~=size(acceptor) | size(acceptor)~=size(fret) )
    error('Data matrix dimensions must agree');
elseif ~isempty(ids) && numel(ids)~=Ntraces
    error('Not enough IDs');
end

if any( isnan(donor(:)) | isnan(acceptor(:)) | isnan(fret(:)) )
    warning('Cannot save NaN values! Converting to zeros');
    donor( isnan(donor(:)) ) = 0;
    acceptor( isnan(acceptor(:)) ) = 0;
    fret( isnan(fret(:)) ) = 0;
end


% Open output file for saving data
fid=fopen(filename,'w');
disp( ['Saving to ' filename] );

% Write time markers (first row) -- universally ignored
fprintf(fid,'%d ', time);
fprintf(fid,'\n');

for j=1:Ntraces
    
    % output name
    name = '';
    if ~isempty(ids)
        name = sprintf('%s ',ids{j});
    end

    % output fluorescence data
    fprintf(fid,'%s',  name);
    fprintf(fid,'%g ', donor(j,:));
    fprintf(fid,'\n');

    fprintf(fid,'%s',  name);
    fprintf(fid,'%g ', acceptor(j,:));
    fprintf(fid,'\n');

    fprintf(fid,'%s',  name);
    fprintf(fid,'%g ', fret(j,:));
    fprintf(fid,'\n');

end % for each trace

fclose(fid);


% END FUNCTION SAVETRACESTXT





function saveTracesQUB( filename, fret )
% FORMAT:
%   All datapoints are concatinated into a M*N column vector;
%   each datapoint is on a new line.
%

fret = fret';

fid=fopen(filename,'w');
disp( ['Saving to ' filename] );

fprintf( fid, '%d\n', fret(:) );

fclose(fid);

% END FUNCTION SAVETRACESQUB





function saveTracesBinary( filename, donor,acceptor,fret, ids, time, metadata )
% FORMAT:
%   struct {
%       int32:
%       int32:   trace length (N)
%       int16:   number of traces (M)
%       string:  IDs (delimited by '-')
%       {int16}: fluorescence data  (M*2 x N matrix)
%       {int32}: time axis (Nx1 vector)
%   }
% 

constants = cascadeConstants;


[nTraces,traceLen] = size(donor);

if ~exist('time','var'),
    time = 1:traceLen; %in frames
end

% Create IDs if not specified
if ~exist('ids','var'),
    [p,name] = fileparts(filename);
    
    ids = cell(Ntraces,1);
    for j=1:Ntraces;
        ids{j} = sprintf('%s_%d', name, j);
    end
end

% Verify input arguments
if any( size(donor)~=size(acceptor) )
    error('Data matrix dimensions must agree');
% elseif exist('ids','var') && numel(ids)~=nTraces
%     error('Not enough IDs');
end

if any( isnan(donor(:)) | isnan(acceptor(:)) | isnan(fret(:)) )
    error('Cannot save NaN values!');
end


% 1) Create IDs if not specified and add to the metadata list
if ~exist('ids','var') || isempty(ids),
    [p,name] = fileparts(filename);
    
    ids = cell(nTraces,1);
    for j=1:nTraces;
        ids{j} = sprintf('%s_%d', name, j);
    end
end

% Remove special characters from IDs
% ids = strrep( ids, '-', '_' );      %- was used as ID seperator, best to avoid.
% ids(' ')='_';      %auto.txt format doesn't allow spaces

metadata(1).ids = sprintf( '%s\t', ids{:} );


% 2) Open file to save data to
fid=fopen(filename,'w');

% 3) Write header data
version = 3; %3 has no data descriptions, 4 does (FUTURE!)
nChannels = 3; %D,A,F

fwrite( fid, 0,         'uint32' );  %identifies the new traces format.
fwrite( fid, 'TRCS',    'char'   );  %format identifier ("magic")
fwrite( fid, version,   'uint16' );  %format version number
fwrite( fid, 9,         'uint8'  );  %trace data type (see loadTraces.m)
fwrite( fid, nChannels, 'uint8'  );

fwrite( fid, [nTraces traceLen], 'uint32' );

% 4) Write fluorescence and FRET data.
fwrite( fid, time,     'single' ); %time axis (in seconds)
fwrite( fid, donor,    'single' );
fwrite( fid, acceptor, 'single' );
fwrite( fid, fret,     'single' );  

% 5) Write metadata pages
fnames = fieldnames( metadata );

for i=1:numel(fnames),
    field = fnames{i};
    metadataText = metadata.(field);
    
    % Write metadata header
    fwrite( fid, numel(field), 'uint8' );
    fwrite( fid, field, 'char' );
    
    fwrite( fid, 0, 'char' ); % 0=char type
    fwrite( fid, numel(metadataText), 'uint32' );
    fwrite( fid, metadataText, 'char' );
end


% Finish up
fclose(fid);





