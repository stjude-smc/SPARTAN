function [donor,acceptor,fret,ids,time,metadata] = loadTraces( ...
                                        filename, constants, indexes )
% LOADTRACES  Loads fluorescence trace files
%
%   [D,A,F,IDs,TIME] = LOADTRACES( FILENAME, CONST, INDEXES )
%   Loads fluorescence donor (D) and acceptor (A) traces;
%   identifiers (ID); and FRET values (F) from the file FILENAME.
%   If the file is .traces (from gettraces), a crosstalk correction is made
%   and FRET is left empty, since there is no such field in these files.
%
%   CONST specifies structure of constants (see cascadeConstants), used
%   only for loading binary files.
%
%   INDEXES specifies the indexes of traces to load. Useful if only a small
%   number of traces are needed from a large file.
%   
%   [D,A,F,IDs,TIME] = LOADTRACES( FILENAME )
%   BINARY FILES ONLY:  Corrections for signal crosstalk and background
%   are made and FRET is calculated.
%

% TODO: make constants the last parameter, since it is rarely used.

% Set empty arguments in case no data is loaded.
[donor,acceptor,fret,time] = deal([]);
ids = {};
metadata = struct();

% If no file is specified, ask for one from the user.
if nargin<1 || isempty(filename),
    [f,p] = uigetfile( {'*.traces','Traces files (*.traces)'; ...
                        '*.txt','Old format traces files (*.txt)'; ...
                        '*.*','All Files (*.*)'}, 'Select a traces file');
    if p==0, return; end
    filename = [p f];
end

if ~exist('indexes','var'), indexes=[]; end
if ~exist('constants','var'), constants = cascadeConstants(); end

% Make sure input file actually exists...
if ~exist(filename,'file'),
    error('File does not exist: %s', filename);
end

[p,n,ext]=fileparts(filename);

if strcmp(ext,'.txt')
    [donor,acceptor,fret,ids,time] = LoadTracesTxt( filename, indexes );
    
elseif strcmp(ext,'.traces')
    [donor,acceptor,fret,ids,time,metadata] = LoadTracesBinary2( ...
                                    filename,constants, indexes );
    
elseif strcmp(ext,'.traces_old')
    [donor,acceptor,fret,ids,time] = LoadTracesBinary_old( ...
                                    filename,constants, indexes );
    
else
    error( ['Filetype not recognized (' ext ')'] );
end

%If IDs are not present, create them
if isempty(ids)
    [p,name] = fileparts(filename);
    name = [name '_'];

    nTraces = size(donor,1);
    ids  = cell(nTraces,1);
    ids(:) = {name};
    strcat( ids, char((1:nTraces)+47)' );

    ids = ids';
end

end %function LoadTraces



%--------------------  LOAD TEXT FILES ------------------------
function [donor,acceptor,fret,ids,time] = LoadTracesTxt( filename, indexes )

% Get file size for calculating how much has been loaded.
d = dir(filename);
fileSize = d.bytes;
clear d;

% Open the traces file
fid=fopen(filename,'r');

% Read the time axis header
time=strread(fgetl(fid),'%f')';
len=length(time);
assert( len>1, 'Cannot parse forQuB file!');


% If the first element of the second line is not a number, then the
% file has molecule ids.
sig=textscan(fid,'%s',1);
sig=sig{:};
hasIDs = isnan(str2double(sig));

% Get the time axis
fseek(fid,-ftell(fid),0);
time=strread(fgetl(fid),'%f')';

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

donor = donor(indexes,:);
acceptor = acceptor(indexes,:);
fret = fret(indexes,:);

% Get trace IDs, if available...
if hasIDs,
    ids = ids(1:3:end-2);
    ids = ids(indexes);
end

% Clean up
close(h);
fclose(fid);
clear Data;

end %function LoadTracesTxt
    
    
function [donor,acceptor,fret,ids,time] = LoadTracesBinary( ...
                                        filename,constants,indexes )

% Open the traces file.
fid=fopen(filename,'r');
len=fread(fid,1,'int32');
Ntraces=fread(fid,1,'int16');

% Read in the trace ids:
c = textscan(fid, '%[^-]', Ntraces/2, 'Delimiter','-');
ids = c{1}';

%If IDs are not present, create them
assert( length(ids) == Ntraces/2, 'LoadTracesBinary: data mismatch' );

% Read in the data:
Data = fread( fid, [Ntraces len], 'int16' );
time = fread( fid,  len, 'int32' );

if isempty(time)
    time = 1:len;
else
    assert( length(time)==len, 'loadTraces: Time axis size mismatch' );
end

% Parse the data into donor, acceptor, etc. arrays.
donor    = double( Data(1:2:end-1,:) );
acceptor = double( Data(2:2:end,  :) );

% Make an adjustment for crosstalk on the camera
acceptor = acceptor - constants.crosstalk*donor;

% Subtract background and calculate FRET
[donor,acceptor,fret] = correctTraces(donor,acceptor,constants,indexes);
ids = ids(indexes);

% Clean up
clear Data;
fclose(fid);

end %function LoadTracesBinary





function [donor,acceptor,fret,ids,time,metadata] = LoadTracesBinary2( ...
                                        filename,constants,indexes )

dataTypes = {'char','uint8','uint16','uint32','uint16', ...
                    'int8', 'int16', 'int32', 'int16', ...
                    'single','double'};  %zero-based
                                    
% 1) Open the traces file.
fid=fopen(filename,'r');

% 2) Read header information
z         = fread( fid, 1, 'uint32' );  %identifies the new traces format.

if z~=0,
    disp('Assuming this is an old format traces file');
    [donor,acceptor,fret,ids,time] = LoadTracesBinary(filename,constants,indexes);
    metadata = struct();
    return;
end

magic     = fread( fid, [1,4], '*char' );  %format identifier ("magic")
version   = fread( fid, 1, '*uint16' );   %format version number
dataType  = fread( fid, 1, '*uint8'  );   %type of trace data - (usually single)
nChannels = fread( fid, 1, '*uint8'  );
nTraces   = fread( fid, 1, 'uint32'  );
traceLen  = fread( fid, 1, 'uint32'  );

% 3) Check validity of header data.
assert( (z==0 && strcmp(magic,'TRCS')), 'loadTraces: invalid header' );
assert( version==3, 'Version not supported!' );

% 4) Read fluorescence and FRET data.
% TODO: To save memory, the data is memory mapped and copied when needed.
assert( nChannels==3, 'Only 2-color FRET data are currently recognized.' );
assert( dataType==9, 'Only single values are supported.');

time = fread( fid, [traceLen,1], 'single' ); %time axis (in seconds)
donor    = fread( fid, [nTraces,traceLen], 'single' );
acceptor = fread( fid, [nTraces,traceLen], 'single' );
fret     = fread( fid, [nTraces,traceLen], 'single' );
szSingle = 4; %size of a single (32 bits)

    
if exist('indexes','var') && ~isempty(indexes),
    donor    = donor(indexes,:);
    acceptor = acceptor(indexes,:);
    fret     = fret(indexes,:);
end

% offset = ftell(fid);
% 
% mmap = memmapfile( filename, 'Offset',offset, 'Format', {...
%                        'single',[nTraces traceLen],'donor'; ...
%                     }, 'Repeat',1,'Writable',false );
% options.forwardToVar='donor';
% donor = mmapPassthrough( mmap, options );
% 
% offset = offset+ (nTraces*traceLen)*szSingle;
% mmap = memmapfile( filename, 'Offset',offset, 'Format', {...
%                        'single',[nTraces traceLen],'acceptor'; ...
%                     }, 'Repeat',1,'Writable',false );
% options.forwardToVar='acceptor';
% acceptor = mmapPassthrough( mmap, options );
% 
% offset = offset+ (nTraces*traceLen)*szSingle;
% mmap = memmapfile( filename, 'Offset',offset, 'Format', {...
%                        'single',[nTraces traceLen],'fret'; ...
%                     }, 'Repeat',1,'Writable',false );
% options.forwardToVar='fret';
% fret = mmapPassthrough( mmap, options );


% 5) Read metadata.
metadata = struct();

while 1,
    % Read the page header.
    titleLength  = fread( fid, 1, 'uint8' );
    title = strtrim(  fread( fid, [1,titleLength], '*char' )  );
    
    pageDatatype = fread( fid, 1, 'uint8' );
    pageSize     = fread( fid, 1, 'uint32' );
    
    %pageY = fread( fid, 1, 'uint32' );
    
    if feof(fid), break; end
    
    % Check validity of field data.
    assert( pageDatatype<numel(dataTypes), 'Invalid field type' );
    if any( isspace(title) ),
        warning('Metadata field titles should not have spaces');
        title(title==' ') = '_';
    end
    
    metadata.(title) = fread( fid, [1,pageSize], ['*' dataTypes{pageDatatype+1}] );
end

% Parse out IDs
if isfield(metadata,'ids')
    c = textscan(metadata.ids, '%s', nTraces, 'Delimiter','\t');
    ids = c{1}';
else
    ids = {};
end


% Clean up
fclose(fid);



end %function LoadTracesBinary






function [donor,acceptor,fret,ids,time] = LoadTracesBinary_old( ...
                                        filename,constants,indexes )
% THIS FUNCTION IS TO LOAD THE DEPRICATED FORMAT.
% Only a few traces files may be found with this format,
% and will have the extension .traces_old...
% These files have no ids!
% NEED TO VERIFY THIS WORKS!
                                    
% Open the traces file.
fid=fopen(filename,'r');
len=fread(fid,1,'int32');
Ntraces=fread(fid,1,'int16');

% Read data
Data = fread( fid, [Ntraces+1 len], 'int16' );
fclose(fid);

% Parse the data into donor, acceptor, etc. arrays.
time     = double( Data(1,:) );
donor    = double( Data(2:2:end-1,:) );
acceptor = double( Data(3:2:end,  :) );
ids = {};

clear Data;

% Make an adjustment for crosstalk on the camera
acceptor = acceptor - constants.crosstalk*donor;

% Subtract background and calculate FRET
[donor,acceptor,fret] = correctTraces(donor,acceptor,constants,indexes);


end %function LoadTracesBinary

