function forOrigin( filename, dwtFilename, outputFilename )
% forOrigin  save trace and idealization data as a text file for easy importing.
%
%   forOrigin() exports trace data in a text format for easy importing into
%   plotting programs like Origin. Traces are listed in order, with 4 columns
%   for each, eg. FRET1, Idl1, Donor1, Acceptor1. The order is designed to make
%   it easy to make "Vertical 2-Panel" plots. The user will be prompted for the
%   input and output files. If no dwell-time file is given, the Idl column will
%   be filled with zeros.
%
%   forOrigin(TRACES, DWT, OUTPUT) uses the supplied file names for the .traces
%   FRET data file, dwell-time file (.dwt), and output .txt file, respectively.
%
%   WARNING: avoiding saving more than 10 traces in this form if possible.
%   Origin may crash with too many columns.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


%% OPTIONS:

% If true, use the average observed FRET value for each trace to determine where
% to draw the idealization line. If false, use the values in the .dwt file.
REESTIMATE = true;


%%

% If no file is specified, ask for one from the user.
if nargin<1,
    filename = getFile('*.traces','Select a traces file');
end
if isempty(filename), return; end  %user hit cancel

% Load traces data
data = loadTraces(filename);
[nTraces,traceLen] = size(data.donor);


% Get idealization filename if not provided.
if nargin<2,
    try
        dwtFilename = findDwt({filename},'raiseError');
        dwtFilename = dwtFilename{1};
    catch
        [p,f] = fileparts(filename);
        dwtFilename = getFile([p f '.qub.dwt'],'Select the corresponding dwell-time file');
    end
end


% Load idealization data (.DWT)
% TODO: what if there are idealizations for multiple FRET traces?
idl = zeros( nTraces,traceLen );
dwells = {};

if exist( dwtFilename, 'file' ),
    % Load idealization from file.
    [dwells,sampling,offsets,model] = loadDWT(dwtFilename);
end

if ~isempty(dwells)
    idl = dwtToIdl( dwells, offsets, traceLen, nTraces );
    
    fretValues = model(:,1);
    nStates = numel(fretValues);
    
    % Convert state sequence to idealization (fret values).
    if REESTIMATE && isChannel(data,'fret'),
        for i=1:size(idl,1),
            % For each state, idealization fret value is average of observed.
            for s=1:nStates,
                fretValues(s) = mean( data.fret(i,idl(i,:)==s) );
            end
            
            % For any states with no dwells, use the global model default
            fretValues(isnan(fretValues)) = model(isnan(fretValues),1);
            
            % Convert idealization from state number to FRET value.
            idl(i, idl(i,:)==0 ) = 1;
            idl(i,:) = fretValues(idl(i,:));
        end
    else
        % Use global FRET values from QuB to create idealization.
        idl( idl==0 ) = 1;
        idl = fretValues(idl);
    end

    % Make sure dimensions are correct for a single trace
    if any( size(idl)==1 ),
        idl = reshape(idl, 1, numel(idl) );
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



% Get field names to save. We want this in a particular order to make it
% easier to make multipanel plots in Origin. Otherwise we could just
% iterate over data.channelNames.
fields = [ dataByName('fret') 'idl' dataByName('donor') dataByName('acceptor') ];
fields = [ fields setdiff(data.channelNames,fields) ]; %get remaining ones


% Store traces data in an output array.
nFields = numel(fields);
totalSize = nTraces*nFields;
output = zeros(traceLen,totalSize+1);

for i=1:nTraces,
    for j=1:numel(fields),
        idx = 1+ (i-1)*nFields + j; %output column index
        if strcmp(fields{j},'idl'),
            output(1:end,idx) = idl(i,:)';
        else
            output(1:end,idx) = data.(fields{j})(i,:)';
        end            
    end
end

output(:,1) = time./1000; %convert to seconds


% Open the output file to save the data.
if nargin<3,
    [p,f] = fileparts(filename);
    outputFilename = fullfile( p, [f '_forOrigin.txt'] );
end

fid = fopen(outputFilename,'w');

% Output header lines
fprintf(fid,'Time (s)');

for i=1:nTraces,
    for j=1:numel(fields),
        fname = fields{j};
        fname(1) = upper(fname(1));
        fprintf( fid, '\t%s_%d', fname,i );
    end
end
fprintf(fid,'\n');

% Output data
for i=1:size(output,1),
    fprintf(fid,'%d\t',output(i,:));
    fprintf(fid,'\n');
end

fclose(fid);






function field = dataByName(fieldname)
% Get all traces with a specified base name (eg, "fret").

field = data.channelNames(  ...
          ~cellfun( @isempty, strfind(data.channelNames,fieldname) )  );
end
      
      
end


