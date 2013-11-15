function forOrigin( filename, dwtFilename, outputFilename )
% Exports trace data in a text format that is easy to import into Origin
% for plotting. Traces are listed in order, with 4 columns for each
% (donor, acceptor, fret, and idealization).
%


% If no file is specified, ask for one from the user.
if nargin<1 || isempty(filename),
    [f,p] = uigetfile( {'*.traces;*.rawtraces','Binary Traces Files (*.traces;*.rawtraces)'; ...
                        '*.txt','Old format traces files (*.txt)'; ...
                        '*.*','All Files (*.*)'}, 'Select a traces file');
    if p==0, return; end
    filename = [p f];
end


% Load traces data
if exist('filename','var'),
    data = loadTraces(filename);
else
    data = loadTraces;
end
[nTraces,traceLen] = size(data.donor);



% Get idealization filename if not provided.
if nargin<2,
    [p,f] = fileparts(filename);
    dwtFilename = [p filesep f '.qub.dwt'];
end

if ~exist(dwtFilename,'file'),
    dwtFilename = '';
    
    [f,p] = uigetfile('*.dwt','Select the corresponding dwell-time file');
    if f~=0,
        dwtFilename = [p f];
    else
        disp('No idealization file found!');
    end
end


% Load idealization data (.DWT)
idl = zeros( nTraces,traceLen );

if ~isempty( dwtFilename ),
    [dwells,sampling,offsets,model] = loadDWT(dwtFilename);

    fretValues = model(:,1);

    idl = dwtToIdl( dwells, traceLen,offsets );
    idl( idl==0 ) = 1;  %assign no-data (no donor) regions to the dark state.
    idl = fretValues(idl);

    % Make sure dimensions are correct for a single trace
    if any( size(idl)==1 ),
        idl = reshape(idl, 1, numel(idl) );
    end

    if size(idl,1)<nTraces,
        idl((size(idl,1)+1):nTraces,:) = 0;
    end
end


% Set time axis if not available,
time = data.time;

if time(1)==1,
    if ~exist('sampling','var')
        f = inputdlg('What is the sampling interval (in ms) for this data?');
        sampling = str2double(f);
    end
    
    tdm = sampling;
    time = 0:tdm:(traceLen*tdm);
    time = time(1:end-1);
end

% Store traces data in an output array
totalSize = nTraces*4;
output = zeros(traceLen,totalSize+1);

for i=1:nTraces,
    idx = 1+ (i-1)*4;
    
    output(1:end,idx+1) = data.fret(i,:)';
    output(1:end,idx+2) = idl(i,:)';
    output(1:end,idx+3) = data.donor(i,:)';
    output(1:end,idx+4) = data.acceptor(i,:)';
end

output(:,1) = time./1000; %convert to seconds


% Output header lines
if nargin<3,
    [p,f] = fileparts(filename);
    outputFilename = [p filesep f '_forOrigin.txt'];
end

fid = fopen(outputFilename,'w');

fprintf(fid,'Time (s)');

for i=1:nTraces,
    fprintf(fid,'\tFRET%d\tIdl%d\tDonor%d\tAcceptor%d',i,i,i,i);
end
fprintf(fid,'\n');

% Output data
for i=1:size(output,1),
    fprintf(fid,'%d\t',output(i,:));
    fprintf(fid,'\n');
end

fclose(fid);
