function data = loadTraces( filename, varargin )
%LOADTRACES  Load fluorescence/FRET data from a .traces format file.
%
%   DATA = LOADTRACES(FILENAME) loads fluorescence/FRET data from a
%   .traces format file into the Traces object DATA. The Traces subclass is
%   chosen based on the channel names.
%
%   DATA = LOADTRACES(..., INDEXES) loads only the specified molecules.
%
%   See also: saveTraces, TracesFret, TracesFret4.

%   FORMAT VERSION HISTORY:
%   1) Simple binary format by JBM starting with the number of traces.
%   2) Extension by JBM to include trace identifier strings.
%   3) First version of this format by DT, including traceMetadata.
%   4) Added channel name list and fileMetadata.
%   5) More flexible metadata format, supporting structs, cell arrays, etc.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


data = [];

dataTypes = {'char','uint8','uint16','uint32','uint64', ...
                    'int8', 'int16', 'int32', 'int64', ...
             'single','double','logical','cell','struct'};  %zero-based
         

% If no file is specified, ask for one from the user.
if nargin<1 || ~ischar(filename),
    filename = getFile;
end
if isempty(filename), return; end


% Load trace data, using distinct functions for each file format.
[~,~,ext]=fileparts(filename);

if strcmp(ext,'.txt')
    data = loadTraces_legacy( filename, varargin{:} );
    return;
    
elseif ~strcmp(ext,'.traces') && ~strcmp(ext,'.rawtraces'),
    warning('Expected .traces or .rawtraces extension.');
end


% Check for obsolete versions of the format.
fid=fopen(filename,'r');

if fread(fid, 1, 'uint32')~=0,  %zero identifies version 3+.
    fclose(fid);
    disp('Assuming this is an old format binary traces file');
    data = loadTraces_legacy(filename, varargin{:});
    return;
end

% Check for the .traces file signature
if ~strcmp(fread(fid, [1,4], '*char'), 'TRCS')
	error('Invalid .traces file header');
end

% Handle old format versions
if fread(fid, 1, '*uint16') < 5,
    fclose(fid);
    data = loadTraces_legacy(filename, varargin{:});
    return;
end


% Read the data dimensions and channel names.
dataType  = fread( fid, 1, '*uint8' );   % data class (9=single)
nChannels = fread( fid, 1, '*uint8' );
nTraces   = fread( fid, 1, 'uint32' );
traceLen  = fread( fid, 1, 'uint32' );

szNames      = fread( fid, 1, 'uint32' );
channelNames = fread( fid, [1,szNames], '*char' );
channelNames = strsplit(channelNames,char(31), 'CollapseDelimiters',false);


% Select subset of traces if indexes are given, silently removing any
% out of range indexes. This is used by traceStat to load data in chunks.
% FIXME: this behavior should not be allowed in the future.
if nargin>=2,
    indexes = varargin{1};
    indexes = indexes( indexes>=1 & indexes<=nTraces );
else
    indexes = 1:nTraces;
end
useAll = isequal(indexes, 1:nTraces);


% Create an empty Traces object.
argin = {numel(indexes), traceLen, channelNames};

if isempty( setdiff(channelNames,{'donor','acceptor','fret'}) ),
    data = TracesFret(argin{:});
    
elseif ismember('donor',channelNames)
    data = TracesFret4(argin{:});
    
elseif ismember('ch1',channelNames),
    data = TracesFluor(argin{:});
    
else
    error('Channel names do not match any available Traces class template');
end


% Read data fields, convert non-float types to double for compatibility.
% Singles are kept as singles to save memory.
if dataType==9,
    star = '*';
else
    star = '';
    disp('Warning: converting non-float data to double!');
end

type = dataTypes{dataType+1};
data.time = fread( fid, [1,traceLen], type );  %time axis (s)

for i=1:nChannels,
    d = fread( fid, [nTraces,traceLen], [star type] );
    
    if useAll,  %faster if indexing is avoided.
        data.(channelNames{i}) = d;
    else
        data.(channelNames{i}) = d(indexes,:);
    end
end


% 5) Read metadata.
root = readMetadata(fid);
fclose(fid);

if isfield(root,'fileMetadata'),
    data.fileMetadata = root.fileMetadata;
end
if isfield(root,'traceMetadata')
    data.traceMetadata = root.traceMetadata(indexes);
end


% Verify all fields are internally consistent and valid.
checkValid(data);


end %function LoadTracesBinary



function metadata = readMetadata(fid)
% M = readMetadata(FID) reads one metadata element from the open file
% handle FID (a .traces file).
%    

% Stars indicate that the data class should be preserved when reading.
dataTypes = {'*char','uint8','uint16','uint32','uint64', ...
                     'int8', 'int16', 'int32', 'int64', ...
             '*single','double','*logical','cell','struct'};

% Read the header, skipping this record if the type is unknown.
dataType = fread( fid, 1, '*uint8'  )+1;   % field type code
nBytes   = fread( fid, 1, '*uint32' );     % content length in bytes
head = ftell(fid);

if dataType<1 || dataType>numel(dataTypes),
    warning('Attempting to skip invalid metadata of type %d.\n', dataType);
    fseek(fid, nBytes, 0);
    metadata = [];
    return;
end


% Read metadata dimensions
ndim       = fread( fid, 1,        '*uint8' );   % number of array dimensions
szMetadata = fread( fid, [1 ndim], '*uint32');   % list of array dimensions
nMetadata  = prod(szMetadata);                   % number of elements.
type = dataTypes{dataType};


% Read metadata contents
if strcmp(type,'cell'),
    metadata = cell( double(szMetadata) );
    for i=1:nMetadata,
        metadata{i} = readMetadata(fid);
    end
    
    
elseif strcmp(type,'struct'),
    metadata = struct();
    nFields = fread(fid, 1, '*uint8');   % number of struct fields
    
    for i=1:nFields
        flen  = fread(fid, 1,       '*uint8');   % field name length
        fname = fread(fid, [1 flen], '*char');   % field name text
        flags = fread(fid, 1,       '*uint8');   % modifiers
        isPacked = bitget(flags,1);
        
        % If not packed, read each struct array element individually.
        if ~isPacked,
            for j=1:nMetadata,
                metadata(j).(fname) = readMetadata(fid);
            end
            continue;
        end
        
        packed = readMetadata(fid);
        if ischar(packed),
            % Unpack strings delimited by char(31).
            packed = strsplit(packed,char(31), 'CollapseDelimiters',false);
            assert( numel(packed)==nMetadata );
            [metadata(1:nMetadata).(fname)] = packed{:};
        
        elseif isnumeric(packed) || islogical(packed),
            % Unpack matrices concatinated along the last dimension.
            if sum(size(packed)>1) == 1,
                temp = num2cell(packed);  %scalar values
            else
                temp = num2cell(packed, 1:ndims(packed)-1);  %matrix values
            end
            [metadata(1:nMetadata).(fname)] = temp{:};
        
        else
            error('Invalid packed datatype %s', class(packed));
        end
    end %for each struct field
    
    
else  %simple numeric, char, logical types
    metadata = fread(fid, nMetadata, type);
end


metadata = reshape(metadata, szMetadata);

% Verify the actual page size matches. (file was corrupted?)
if ftell(fid)-head ~= nBytes,
    error('Metadata page size mismatch');
end


end %FUNCTION readMetadata






