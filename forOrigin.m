function forOrigin( filename, dwtFilename )
% TODO: only assign low-FRET state as zero if feature is enabled;
% and only do so for a maximum of X frames.

maxBlink = 5; %frames


% Load traces data
if exist('filename','var'),
    data = loadTraces(filename);
else
    data = loadTraces;
end
[nTraces,traceLen] = size(data.donor);

% Load idealization data (.DWT)
idl = zeros( nTraces,traceLen );

if ~exist('dwtFilename','var'),
    [f,p] = uigetfile('*.dwt');
    if f~=0,
        dwtFilename = [p f];
        [dwells,sampling,offsets,model] = loadDWT(dwtFilename);
        
        fretValues = model(:,1);
        %fretValues(fretValues==0.1) = 0.01; %for LeuT?

        idl = dwtToIdl( dwells, traceLen,offsets );
        
        idl = idl';
        blinks = rleFilter( idl==1, maxBlink );
        idl( idl==1 ) = 2; %ignore blinks
        idl( idl==0 ) = 1; 
        idl( blinks ) = 1; %restore long blinks
        idl = idl';
        
        idl = fretValues(idl);
        
        % Make sure dimensions are correct for a single trace
        if any( size(idl)==1 ),
            idl = reshape(idl, 1, numel(idl) );
        end
        
        if size(idl,1)<nTraces,
            idl((size(idl,1)+1):nTraces,:) = 0;
        end
    else
        disp('No idealization file found!');
    end
end


% Set time axis if not available,
time = data.time;

if time(1)==1,
    if ~exist('sampling','var')
        f = inputdlg('What is the sampling interval (in ms) for this data?');
        sampling = str2double(f)
    end
    
    tdm = sampling;
    time = 0:tdm:(traceLen*tdm);
    time = time(1:end-1);
else
    sampling = time(2)-time(1);
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
fid = fopen('traces.txt','w');

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
