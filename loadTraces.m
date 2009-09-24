function [donor,acceptor,fret,ids,time] = loadTraces( ...
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

if nargin<1,
    [f,p] = uigetfile('*.txt','Select a traces file');
    if f==0, return; end
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
    [donor,acceptor,fret,ids,time] = LoadTracesBinary( ...
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

% Clean up
clear Data;
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

