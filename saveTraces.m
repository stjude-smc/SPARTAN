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
%   {
%       uint32:  zero (distinguishes from older format without a header)
%       char4:   TRCS (magic string to unambiguously identify file)
%       unit16:  version number (4)
%       unit8:   traces data type (9=single, zero-based)
%       
%       uint8:   number of channels (C)
%       uint32:  number of traces (M)
%       uint32:  number of frames per traces (N)
%
%       unit32:  length of channel names string...
%       {char}:  list of channel names (delimited by char/31)
%
%       {single}: fluorescence/fret data  (C x M x N matrix)
%       {uint32}: time axis (Nx1 vector)
%
%       For each metadata page (until end-of-file):
%         uint32:  field title length
%         {char}:  filed title
%         uint8:   data type identifider (see dataTypes)
%         uint32:  metadata length (in units, not bytes)
%         {xxx}:   metadata content
%   }
% 
% There are two groups of metadata pages: fileMetadata and traceMetadata.
% fileMetadata are arbitrary-length data that applies to all traces in the
% file (e.g., power meter reading, EM gain settings, etc). traceMetadata is
% a structure array with one element per trace that contains data specific
% to each trace (e.g., molecule locations in the field-of-view, molecule
% identifiers, etc). These sections are delineated by special
% "section_heading" metadata pages containing the section name
% (e.g., traceMetadata).
% 

global dataTypes;
dataTypes = {'char','uint8','uint16','uint32','uint16', ...
                    'int8', 'int16', 'int32', 'int16', ...
                    'single','double'};  %zero-based

[nTraces,traceLen] = size(data.donor);

if ~isfield(data,'time'),
    data.time = 1:traceLen; %in frames
end


% Legacy code (or laziness) support. If no channel names are given, try to
% guess based on what data are present.
if ~isfield(data,'channelNames'),
    data.channelNames = {'donor','acceptor','fret'};
end


% Verify input arguments
nChannels = numel(data.channelNames);

for i=1:nChannels,
    ch = data.( data.channelNames{i} );
    if any( size(ch) ~= size(data.donor) ),
        error('Data matrix dimensions must agree');
    end
    if any( isnan(ch) )
        error('Cannot save NaN values!');
    end
end


[p,f,e] = fileparts(filename);
assert( ~isempty(strfind(e,'traces')), 'Binary format traces files must have a ".*traces" extension' );


% 1) Create IDs if not specified and add to the metadata list
if ~isfield(data,'traceMetadata');
    data.traceMetadata = struct();
end

if ~isfield(data.traceMetadata,'ids'),
    disp('saveTraces: no ids detected. Recreating them.');
    for i=1:nTraces;
        data.traceMetadata(i).ids = sprintf('%s#%d', filename, i);
    end
end


% 2) Open file to save data to
fid=fopen(filename,'w');

% 3) Write header data
version = 4; % ver 4 adds file-global metadata (see #6 below).

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
    
    % Write section header to identify this as per-trace metadata.
    writeMetadata( fid, 'section_heading', 'traceMetadata' );
end

for i=1:numel(fnames),
    fname = fnames{i};
    m = data.traceMetadata(1).(fname); %just get first element for determining data type
    
    % Collapse strucutre array into a single field for serialization.
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
    
    writeMetadata( fid, fname, metadata );
end
 

% 5) Write global metadata pages (if any).
% Note that here these are not arrays over traces (the structure is not
% an array), so it is simpler to process than traceMetadata.
fnames = {};
if isfield(data,'fileMetadata') && numel( data.fileMetadata )>0,
    fnames = fieldnames( data.fileMetadata );
    
    % Write section header to identify this as per-trace metadata.
    writeMetadata( fid, 'section_heading', 'fileMetadata' );
end

for i=1:numel(fnames),
    fname = fnames{i};
    m = data.fileMetadata.(fname);
    
    % Collapse strucutre array into a single field for serialization.
    if ~isnumeric(m) || ischar(m),
        warning( 'saveTraces:badMetadataType', ['Unsupported metadata field type: ' fname ' (' class(m) ')'] );
    else
        writeMetadata( fid, fname, m );
    end
end


% Finish up
fclose(fid);




function writeMetadata( fid, fname, metadata )
% Writes a single metadata "page" to currently open traces file
% 
% TODO(?): allow cell arrays of strings to be written as delimited lists.


% Determine data type code and verify it is an allowed type.
% A warning is ok here since this metadata info will just be skipped.
global dataTypes;
fieldDataType = find( strcmp(class(metadata),dataTypes) );
if isempty( fieldDataType ),
    warning( 'saveTraces:badMetadataType', ['Unsupported metadata field type: ' fname ' (' class(metadata) ')'] );
    return;
end

% Write metadata header and content.
fwrite( fid, numel(fname), 'uint8' );       % field title length
fwrite( fid, fname, 'char' );               % field title text
fwrite( fid, fieldDataType-1, 'uint8' );    % field type code
fwrite( fid, numel(metadata), 'uint32' );   % content length
fwrite( fid, metadata, class(metadata) );   % content


% end function writeMetadata




