function saveTraces( filename, varargin )
% SAVETRACES    Saves trace data to file
%
%    SAVETRACES( FNAME, [FORMAT,] DATA )
%
% Saves fluorescence and FRET data to file. FNAME is the filename to save
% to, FORMAT is the file format, and data are is a structure with all of
% the fluorescence/FRET traces. 'txt' format saves the time axis, IDs,
% donor/acceptor/fret traces as plain text. 'qub' format saves just the
% FRET data as text (for importing into QuB software). 'traces' format
% saves as the binary traces format (default). If format is not given, the
% standard 'traces' format is assumed.
%
% DATA typically contains the following fields: (see gettraces.m)
%   - channelNames (cell array of strings: donor, acceptor fret, ...)
%   - donor (donor fluorescence)
%   - acceptor (acceptor fluorescence)
%   - fret (fret ratio as A/(A+D))
%   - factor (miscellaneous fluorescence signal, eg, factor binding, optional)
%   - traceMetadata (structure array of metadata for each molecule, optional)
%   - fileMetadata (structure of metadata that is applies to the whole file, optional)
%
% The channel names can be anything really, but these are values that most
% of the code expects. Multi-pair FRET not supported yet. FIXME
% 
% For 'qub' files, data may just be the FRET traces
% 

% If no format name given, assume traces file.
if nargin==2 && isstruct( varargin{1} ),
    saveTracesBinary( filename, varargin{1} );
    return;
end

assert( nargin==3, 'Invalid number of arguments' );

% Otherwise, determine the data type by the format.
format = varargin{1};
data   = varargin{2};

switch format
    case 'txt'
        saveTracesTxt( filename, data );
        
    case 'traces'
        saveTracesBinary( filename, data );
        
    case 'qub'
        if isstruct(data)
            fret = data.fret;
        else
            fret = data;
        end
        
        saveTracesQUB( filename, fret );
        
    otherwise
        error('Unknown file format');
end

% END FUNCTION SAVETRACES



function saveTracesTxt( filename, data )
% FORMAT:
%   1 2 3 4 ... N
%   ID1 donor data
%   ID1 acceptor data
%   ID1 FRET data
%   ID2 donor data
%   ...
% 

[p,f,e] = fileparts(filename);
assert( ~isempty(strcmp(e,'.txt')), 'TXT format traces files must have a ".txt" extension' );

assert( isfield(data,'donor') & isfield(data,'acceptor') & isfield(data,'fret'), ...
        'Data to save must include, donor, acceptor, and fret traces' );

[Ntraces,tlen] = size(data.donor);

if ~isfield(data,'time'),
    time = 1:tlen; %in frames
end

% Create IDs if not specified
if ~isfield(data,'ids'),
    [p,name] = fileparts(filename);
    
    ids = cell(Ntraces,1);
    for j=1:Ntraces;
        ids{j} = sprintf('%s_%d', name, j);
    end
else
    ids = data.ids;
end

% Remove special characters from IDs
ids = strrep( ids, '-', '_' );      %- is used as ID seperator, must be avoided
ids = strrep( ids, ' ', '_' );      %auto.txt format doesn't allow spaces

% Verify input arguments
if any( size(data.donor)~=size(data.acceptor) | size(data.acceptor)~=size(data.fret) )
    error('Data matrix dimensions must agree');
elseif ~isempty(ids) && numel(ids)~=Ntraces
    error('Not enough IDs');
end

if any( isnan(data.donor(:)) | isnan(data.acceptor(:)) | isnan(data.fret(:)) )
    warning('Cannot save NaN values! Converting to zeros');
    data.donor( isnan(data.donor(:)) ) = 0;
    data.acceptor( isnan(data.acceptor(:)) ) = 0;
    data.fret( isnan(data.fret(:)) ) = 0;
end


% Open output file for saving data
fid=fopen(filename,'w');
disp( ['Saving to ' filename] );

% Write time markers (first row) -- universally ignored
fprintf(fid,'%d ', data.time);
fprintf(fid,'\n');

for j=1:Ntraces
    
    % output name
    name = '';
    if ~isempty(ids)
        name = sprintf('%s ',ids{j});
    end

    % output fluorescence data
    fprintf(fid,'%s',  name);
    fprintf(fid,'%g ', data.donor(j,:));
    fprintf(fid,'\n');

    fprintf(fid,'%s',  name);
    fprintf(fid,'%g ', data.acceptor(j,:));
    fprintf(fid,'\n');

    fprintf(fid,'%s',  name);
    fprintf(fid,'%g ', data.fret(j,:));
    fprintf(fid,'\n');

end % for each trace

fclose(fid);


% END FUNCTION SAVETRACESTXT





function saveTracesQUB( filename, fret )
% FORMAT:
%   All datapoints are concatinated into a M*N column vector;
%   each datapoint is on a new line.
%

[p,f,e] = fileparts(filename);
assert( ~isempty(strcmp(e,'.txt')), 'QuB TXT format traces files must have a ".txt" extension' );

fret = fret';

fid=fopen(filename,'w');
disp( ['Saving to ' filename] );

fprintf( fid, '%d\n', fret(:) );

fclose(fid);

% END FUNCTION SAVETRACESQUB





function saveTracesBinary( filename, data )
% Saves trace data in binary format. Required fields: donor, acceptor, fret
%
% FORMAT:
%   struct {
%       (header data -- see below)
%       uint8:   number of channels (C)
%       uint32:  number of traces (M)
%       uint32:  number of frames per traces (N)
%
%       {single}: fluorescence/fret data  (C x M x N matrix)
%       {int32}: time axis (Nx1 vector)
%
%       For each metadatafield (until end-of-file):
%         uint32:  field title length
%         {char}:  filed title
%         uint8:   data type identifider (see dataTypes)
%         uint32:  metadata length (in units, not bytes)
%         {xxx}:   metadata content
%   }
% 

assert( isfield(data,'donor') & isfield(data,'acceptor') & isfield(data,'fret'), ...
        'Data to save must include, donor, acceptor, and fret traces' );

dataTypes = {'char','uint8','uint16','uint32','uint16', ...
                    'int8', 'int16', 'int32', 'int16', ...
                    'single','double'};  %zero-based

[nTraces,traceLen] = size(data.donor);

if ~isfield(data,'time'),
    data.time = 1:traceLen; %in frames
end


% Verify input arguments
if any( size(data.donor)~=size(data.acceptor) )
    error('Data matrix dimensions must agree');
end

if any( isnan(data.donor(:)) | isnan(data.acceptor(:)) | isnan(data.fret(:)) )
    error('Cannot save NaN values!');
end

[p,f,e] = fileparts(filename);
assert( ~isempty(strfind(e,'traces')), 'Binary format traces files must have a ".*traces" extension' );


% 1) Create IDs if not specified and add to the metadata list
if ~isfield(data,'traceMetadata');
    data.traceMetadata = struct();
end

if ~isfield(data.traceMetadata,'ids'),
    for i=1:nTraces;
        data.traceMetadata(i).ids = sprintf('%s#%d', filename, i);
    end
end


% 2) Open file to save data to
fid=fopen(filename,'w');

% 3) Write header data
version = 4; % ver 4 adds file-global metadata (see #6 below).
nChannels = numel(data.channelNames);

fwrite( fid, 0,         'uint32' );  %identifies the new traces format.
fwrite( fid, 'TRCS',    'char'   );  %format identifier ("magic")
fwrite( fid, version,   'uint16' );  %format version number
fwrite( fid, 9,         'uint8'  );  %trace data type (single, see loadTraces.m)
fwrite( fid, nChannels, 'uint8'  );
fwrite( fid, [nTraces traceLen], 'uint32' );

% Write channel names (donor, acceptor, fret, etc).
channelNames = strcat( data.channelNames, char(31) );
channelNames = strcat( channelNames{:} );
channelNames = channelNames(1:end-1); %removing trailing seperator

fwrite( fid, numel(channelNames), 'uint32' );
fwrite( fid, channelNames, 'char' );


% 4) Write fluorescence and FRET data.
fwrite( fid, data.time, 'single' ); %time axis (in seconds)

for i=1:nChannels,
    fwrite( fid, data.(data.channelNames{i}), 'single' );
end


% 5) Write per-trace metadata pages (if any)
fnames = {};
if isfield(data,'traceMetadata') && numel( data.traceMetadata )>0,
    fnames = fieldnames( data.traceMetadata );
end

for i=1:numel(fnames),
    fname = fnames{i};
    
    % Collapse strucutre array into a single field for serialization.
    m = data.traceMetadata(1).(fname);
    
    if isnumeric(m),
        metadata = [data.traceMetadata.(fname)];
    elseif ischar(m),
        metadata = strcat( {data.traceMetadata.(fname)}, char(31) );
        metadata = [ metadata{:} ];
        metadata = metadata(1:end-1); %removing trailing seperator
    else
        warning( 'saveTraces:badMetadataType', ['Unsupported metadata field type: ' fname ' (' class(m) ')'] );
        continue;
    end
    
    % Write metadata header
    fwrite( fid, numel(fname), 'uint8' );  % field title length
    fwrite( fid, fname, 'char' );          % field title text
    
    fieldDataType = find( strcmp(class(metadata),dataTypes) );
    if isempty( fieldDataType ),
       error( 'Unsupported metadata field data type' ); 
    end
    
    % Write field content
    fwrite( fid, fieldDataType-1, 'uint8' );
    fwrite( fid, numel(metadata), 'uint32' );
    fwrite( fid, metadata, class(metadata) );
end
 

% 5) Write global metadata pages (if any).
fnames = {};
if isfield(data,'fileMetadata') && numel( data.fileMetadata )>0,
    fnames = fieldnames( data.fileMetadata );
end

for i=1:numel(fnames),
    % Get next field.
    % Note that here these are not arrays over traces (the structure is not
    % an array), so it is simpler to process than traceMetadata.
    fname = fnames{i};
    m = data.fileMetadata.(fname);
    
    if ~isnumeric(m) && ~ischar(m),
        warning( 'saveTraces:badMetadataType', ['Unsupported metadata field type: ' fname ' (' class(m) ')'] );
        continue;
    end
    
    % Write metadata header
    fwrite( fid, numel(fname), 'uint8' );  % field title length
    fwrite( fid, fname, 'char' );          % field title text
    
    fieldDataType = find( strcmp(class(metadata),dataTypes) );
    if isempty( fieldDataType ),
       error( 'saveTraces:badMetadataType', ['Unsupported metadata field type: ' fname ' (' class(m) ')'] );
    end
    
    % Write field content
    fwrite( fid, fieldDataType-1, 'uint8' );
    fwrite( fid, numel(metadata), 'uint32' );
    fwrite( fid, metadata, class(metadata) );
end


% Finish up
fclose(fid);





