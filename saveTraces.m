function saveTraces( filename, format, varargin )
% SAVETRACES    Saves trace data to file
%
%    SAVETRACES( FNAME, 'txt',    D,A,F,IDs, time )
%    SAVETRACES( FNAME, 'traces', D,A,  IDs, time )
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
    error('Cannot save NaN values!');
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





function saveTracesBinary( filename, donor,acceptor, ids, time )
% FORMAT:
%   struct {
%       int32:   trace length (N)
%       int16:   number of traces (M)
%       string:  IDs (delimited by '-')
%       {int16}: fluorescence data  (M*2 x N matrix)
%       {int32}: time axis (Nx1 vector)
%   }
% 

[Ntraces,tlen] = size(donor);

if ~exist('time','var'),
    time = 1:tlen; %in frames
end

% Verify input arguments
if any( size(donor)~=size(acceptor) )
    error('Data matrix dimensions must agree');
% elseif ~exist('ids','var') && numel(ids)~=Ntraces
%     error('Not enough IDs');
end

if any( isnan(donor(:)) | isnan(acceptor(:)) )
    error('Cannot save NaN values!');
end


% Create IDs if not specified
if ~exist('ids','var') || isempty(ids),
    [p,name] = fileparts(filename);
    
    ids = cell(Ntraces,1);
    for j=1:Ntraces;
        ids{j} = sprintf('%s_%d', name, j);
    end
end

% Remove special characters from IDs
ids = strrep( ids, '-', '_' );      %- is used as ID seperator, must be avoided
ids = strrep( ids, ' ', '_' );      %auto.txt format doesn't allow spaces


% Open file to save data to
fid=fopen(filename,'w');

% Save header data
fwrite(fid,tlen,   'int32');
fwrite(fid,Ntraces*2,'int16');

for j=1:Ntraces
    fprintf(fid, '%s-', ids{j});
end

% Save trace data
traces = zeros(Ntraces*2,tlen);
traces(1:2:end,:) = donor;
traces(2:2:end,:) = acceptor;

fwrite(fid,traces,'int16');
fwrite(fid,time,'int32');

% Finish up
fclose(fid);


% fid=fopen(save_fname,'w');
% 
% fwrite(fid,tlen,'int32');
% fwrite(fid,Ntraces,'int16');  % max of 32,000 traces!
% 
% name=strrep(name,'-','_');  %- is used as ID seperator, must be avoided
% name=strrep(name,' ','_');  %auto.txt format doesn't allow spaces
% 
% for j=1:Ntraces/2;
%     fprintf(fid, '%s_%d-', name, j);
% end
% fwrite(fid,traces,'int16');
% fclose(fid);
% clear traces;













