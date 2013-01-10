function data = loadTraces( filename, indexes )
% LOADTRACES  Loads fluorescence trace files
%
%   DATA = LOADTRACES( FILENAME )  loads fluorescence DATA and metadata
%   from file. For non-obselete versions of this format, no corrections
%   are made here (they should have been done in gettraces.m).
%
%   DATA = LOADTRACES( FILENAME, INDEXES )  loads only specified traces.
%
%   For two-color FRET experiments, fields include:
%      donor, acceptor, fret, time
%   Where donor and acceptor are fluorescence data, fret is calculated FRET
%   ratio trace, time is the time axis (in seconds).
% 
%   Metadata fields include:
%      channelNames:   names of the data channels ('donor','acceptor',etc)
%      traceMetadata:  data specific to each trace (structure array)
%      fileMetadata:   data specific to each the entire file
%   
%   Additional fields may be present for addition fluorescence channels.
%   This includes 'factor' for factor binding or 'acceptor2' for a
%   three-color FRET experiment with a single donor and two acceptors.
%   Four-color FRET and more complex experiments are not supported. FIXME
%

data = struct();

% If no file is specified, ask for one from the user.
if nargin<1 || isempty(filename),
    [f,p] = uigetfile( {'*.traces;*.rawtraces','Binary Traces Files (*.traces;*.rawtraces)'; ...
                        '*.txt','Old format traces files (*.txt)'; ...
                        '*.*','All Files (*.*)'}, 'Select a traces file');
    if p==0, return; end
    filename = [p f];
end

if ~exist('indexes','var'), indexes=[]; end

% Make sure input file actually exists...
if ~exist(filename,'file'),
    error('File does not exist: %s', filename);
end

% Load trace data, using distinct functions for each file format.
[p,n,ext]=fileparts(filename);

if strcmp(ext,'.txt')
    data = LoadTracesTxt( filename, indexes );
    
elseif strcmp(ext,'.traces') || strcmp(ext,'.rawtraces')
    data = LoadTracesBinary2( filename, indexes );
    
else
    error( ['Filetype not recognized (' ext ')'] );
end

% If IDs are not present, create them.
if ~isfield(data,'traceMetadata') || ~isfield(data.traceMetadata,'ids')
    for i=1:size(data.donor,1),
        data.traceMetadata(i).ids = sprintf( '%s#%d', filename, i );
    end
end

end %function LoadTraces



%--------------------  LOAD TEXT FILES ------------------------
function data = LoadTracesTxt( filename, indexes )

% Get file size for calculating how much has been loaded.
d = dir(filename);
fileSize = d.bytes;
clear d;

% Open the traces file
fid=fopen(filename,'r');

% Read the time axis header
data.time=strread(fgetl(fid),'%f')';
len=length(data.time);
assert( len>1, 'Cannot parse forQuB file!');


% If the first element of the second line is not a number, then the
% file has molecule ids.
sig=textscan(fid,'%s',1);
sig=sig{:};
hasIDs = isnan(str2double(sig));

% Get the time axis
fseek(fid,-ftell(fid),0);
data.time=strread(fgetl(fid),'%f')';

% Extract intensity information (and IDs) from file.
h = waitbar(0,'Loading trace data...');

ids = cell(0,1);
i = 1;
Data = cell(0);
while 1
    if hasIDs
        id = textscan(fid, '%s',1);
        if isempty(id{1}), break; end  %end of file

        ids{i} = id{1}{1};
    end
    
    line = textscan(fid, '%f', len);
    if isempty(line{1}), break; end  %end of file
    
    Data{i} = line{1};
    i = i+1;
    
    % update wait bar occasionally based on amount of data read.
    if mod(i,20)==0,
        waitbar( 0.99*ftell(fid)/fileSize, h );
    end
end
Data = cell2mat(Data)';

% Split data into donor, acceptor, and FRET channels.
donor=Data(1:3:end-2,:);
acceptor=Data(2:3:end-1,:);
fret  = Data(3:3:end,:);
assert( all( size(fret)==size(donor)) );

% Select only molecules specified
if nargin<2 || isempty(indexes),
    indexes = 1:size(donor,1);
end

data.donor = donor(indexes,:);
data.acceptor = acceptor(indexes,:);
data.fret = fret(indexes,:);

% Get trace IDs, if available...
if hasIDs,
    ids = ids(1:3:end-2);
    ids = ids(indexes);
    data.traceMetadata = struct( 'ids', ids );
end

% Clean up
close(h);
fclose(fid);
clear Data;

end %function LoadTracesTxt
    
    
function data = LoadTracesBinary( filename, indexes )
% This loads binary trace data from the old format, which did not include
% metadata (obsolete as of 8/15/2012).

constants = cascadeConstants();

% Open the traces file.
fid=fopen(filename,'r');
len=fread(fid,1,'int32');
Ntraces=fread(fid,1,'int16');

% Read in the trace ids (required!)
c = textscan(fid, '%[^-]', Ntraces/2, 'Delimiter','-');
ids = c{1}';
assert( length(ids) == Ntraces/2, 'LoadTracesBinary: data mismatch' );

% Read in the data:
Input = fread( fid, [Ntraces len], 'int16' );
data.time = fread( fid,  len, 'int32' );

if isempty(data.time)
    data.time = 1:len;
else
    assert( length(data.time)==len, 'loadTraces: Time axis size mismatch' );
end

% Parse the data into donor, acceptor, etc. arrays.
donor    = double( Input(1:2:end-1,:) );
acceptor = double( Input(2:2:end,  :) );

% Make an adjustment for crosstalk on the camera
acceptor = acceptor - constants.crosstalk*donor;

% Subtract background and calculate FRET
[data.donor,data.acceptor,data.fret] = correctTraces( ...
                                donor,acceptor,constants,indexes);
data.traceMetadata = struct( 'ids', ids(indexes) );

% Clean up
clear Data;
fclose(fid);


end %function LoadTracesBinary





function data = LoadTracesBinary2( filename, indexes )
% Loads binary data from the current standard trace format.

dataTypes = {'char','uint8','uint16','uint32','uint16', ...
                    'int8', 'int16', 'int32', 'int16', ...
                    'single','double'};  %zero-based
                                    
% 1) Open the traces file.
fid=fopen(filename,'r');

% 2) Read header information
z = fread( fid, 1, 'uint32' );  %identifies the new traces format.

if z~=0,
    disp('Assuming this is an old format binary traces file');
    data = LoadTracesBinary(filename,indexes);
    return;
end

magic     = fread( fid, [1,4], '*char' );  %format identifier ("magic")
version   = fread( fid, 1, '*uint16' );   %format version number
dataType  = fread( fid, 1, '*uint8'  );   %type of trace data - (usually single)
nChannels = fread( fid, 1, '*uint8'  );
nTraces   = fread( fid, 1, 'uint32'  );
traceLen  = fread( fid, 1, 'uint32'  );

data.nChannels = nChannels;

% 3) Check validity of header data.
assert( (z==0 && strcmp(magic,'TRCS')), 'loadTraces: invalid header' );
assert( version>=3 && version<=4, 'Version not supported!' );


% 4a) Read data field names
if version>3, %field names are specified in the new format.
    szNames = fread( fid, 1, 'uint32' );
    channelNames = strtrim(  fread( fid, [1,szNames], '*char' )  );

    c = textscan(channelNames, '%s', 'Delimiter',char(31));
    channelNames = c{1}';
else
    % For backward compatibility (with ver. 3), assume 2-color FRET.
    assert( nChannels==3 );
    channelNames = {'donor','acceptor','fret'};
end

data.channelNames = channelNames;


% 4b) Read fluorescence and FRET data.
% TODO: To save memory, the data is memory mapped and copied when needed.
assert( dataType==9, 'Only single values are supported.');

data.time = fread( fid, [traceLen,1], 'single' ); %time axis (in seconds)

% Select subset of traces if indexes are given.
if nargin<2 || isempty(indexes),
    indexes = 1:nTraces;
end

% Read data fields.
for i=1:nChannels,
    d = fread( fid, [nTraces,traceLen], 'single' );
    data.( channelNames{i} ) = d(indexes,:);
end


% 5) Read metadata.
% Adding this just to give the traceMetadata struct a basic structure,
% which is needed for building structure arrays. 'temp' is removed later.
traceMetadata = struct( 'temp', num2cell(indexes) );
mode = 1; %1=reading traceMetadata (default), 2=reading fileMetadata

while 1,
    % Read the page header.
    titleLength  = fread( fid, 1, 'uint8' );
    title = strtrim(  fread( fid, [1,titleLength], '*char' )  );
    
    pageDatatype = fread( fid, 1, 'uint8' );
    pageSize     = fread( fid, 1, 'uint32' );
    
    if feof(fid), break; end
    
    % Check validity of field data.
    assert( pageDatatype<numel(dataTypes), 'Invalid field type' );
    if any( isspace(title) ),
        warning('Metadata field titles should not have spaces');
        title(title==' ') = '_';
    end
    
    m = fread( fid, [1,pageSize], ['*' dataTypes{pageDatatype+1}] );
    
    % Extract ids (delimited text).
    % FIXME: do this for any field with char(31) -- delimited text.
    if strcmp(title,'ids'),
        c = textscan(m, '%s', nTraces, 'Delimiter',char(31));
        m = c{1}';
    end
    
    % Look for section markers that delimit traceMetadata from fileMetadata
    if strcmp(title,'section_heading'),
        assert( ischar(m) );
        if strcmp(m,'traceMetadata'),
            mode=1;
        elseif strcmp(m,'fileMetadata'),
            mode=2;
        else
            error('Unknown section seperator');
        end
    end
    
    % Convert into structure array for traceMetadata.
    if mode==1,
        if iscell(m),
            [traceMetadata.(title)] = deal( m{indexes} );
        elseif isnumeric(m),
            d = num2cell(m(indexes));
            [traceMetadata.(title)] = deal( d{:} );
        end
    
    % Simply add fields (which may be arrays) for fileMetadata.
    elseif mode==2,
        data.fileMetadata.(title) = m;
    end
end


% Clean up
data.traceMetadata = rmfield(traceMetadata,'temp');
fclose(fid);


end %function LoadTracesBinary

