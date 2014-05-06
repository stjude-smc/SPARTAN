function output = combineDatasets( filenames, outFilename )
% combineDatasets  Combine several smFRET data files into one file.
%
%  combineDatasets( FILES, OUTPUT )
%  Each smFRET data file (see loadTraces.m) specified in the the cell array
%  FILES is loaded and combined into a single, large dataset. The resulting
%  dataset is saved to the OUTPUT filename. If files are of differing
%  lengths, they are truncated to the minimal size. If input arguments are
%  not specified, the user will be prompted for them.
%


%% Get file names if not specified.
if nargin<1,
    filenames = getFiles;
    if isempty(filenames), return; end;
end

if nargin<2,
    [f,p] = uiputfile('*.traces','Select output filename','combined.traces');
    if f==0, return; end
    outFilename = [p f];
end

nFiles = numel(filenames);
if nFiles<1, return; end


h = waitbar(0,'Combining datasets');


%% Load data
nTraces  = zeros(1,nFiles);
traceLen = zeros(1,nFiles);

traceMetadataFields = {}; %field names that are common across files.
channelNames = {}; %channel names that are common across all files.

dwt=cell(nFiles,1);
sampling=zeros(nFiles,1);
offsets=cell(nFiles,1);
model=cell(nFiles,1);

for i=1:nFiles,

    % Load traces from file. If any of the files have a slightly different
    % structure (fields missing or in a different order), this will fail!
    d = loadTraces(filenames{i});
    [nTraces(i),traceLen(i)] = size(d.donor);
    assert( traceLen(i)>1 );
    
    % Select only the fields that are common to both files. Generally
    % all fields will be in common. This is just a precaution.
    if i>1,

        fieldsToRemove = setxor(fieldnames(d),fieldnames(data));
        
        if ~isempty(fieldsToRemove),
            d    = rmfield( d,    setdiff(fieldnames(d),fieldnames(data)) );
            data = rmfield( data, setdiff(fieldnames(data),fieldnames(d)) );
            
            text = sprintf( '%s, ',fieldsToRemove{:} );
            warning( 'combineDatasets:missingDataFields', ...
                     'Some fields were removed because they are not in all files: %s', ...
                     text(1:end-2)  );
        end
    end
    data(i) = d;
    
    % Get common field names for data and metadata fields.
    if i==1,
        channelNames = d.channelNames;
        traceMetadataFields = fieldnames(d.traceMetadata);
    else
        % FIXME: channels may be in a different order. This will give an
        % error here, but it is possible to reorder the data.
        assert( all(strcmp(channelNames,d.channelNames)), ...
                                      'Data channels must be the same' );
        %if ~isempty( setxor(channelNames,d.channelNames) ),
        %    close(h);
        %    error('All files must have the same data channels!');
        %end
        
        %channelNames = intersect(channelNames,d.channelNames);
        traceMetadataFields = intersect( ...
                        traceMetadataFields, fieldnames(d.traceMetadata) );
    end
    
    
    % Verify the data are valid. We can't save NaN values to file!
    for c=1:numel(channelNames),
        chName = channelNames{c};
        assert(  ~any( isnan(data(i).(chName)(:)) )  );
    end
    
    
    % Look for an idealization. These will be combined if every dataset has one.
    [path,file] = fileparts( filenames{i} );
    dwt_fname = [path filesep file '.qub.dwt'];
    if exist(dwt_fname,'file'), 
        [dwt{i},sampling(i),offsets{i},model{i}] = loadDWT( dwt_fname );
    end
    
    waitbar(0.7*i/nFiles,h);
end


% Get the longest time axis
[~,longest] = max(traceLen);
time = data(longest).time;


%% Resize traces so they are all the same length
newTraceLength = traceLen(1);

if min(traceLen) ~= max(traceLen),
    
    choice = questdlg('Traces are different lengths. How do you want to modify the traces so they are all the same length?', ...
                  'Resize traces', 'Truncate','Extend','Cancel', 'Cancel');

    % ---- User hit cancel, do nothing.
    if strcmp(choice,'Cancel'),
        close(h);
        return;
        
    % ---- Truncate traces to the same length
    elseif strcmp(choice,'Truncate'),
        
        newTraceLength = min(traceLen);

        for i=1:nFiles,
            % Resize each channel in this file.
            for c=1:numel(channelNames),
                chName = channelNames{c};
                data.(chName) = data.(chName)(:,1:newTraceLength);
            end
            
            % Convert dwell-times to an idealization, which are easy to
            % truncate, and convert back for saving later.
            if ~isempty(dwt{i}),
                idl = dwtToIdl( dwt{i}, traceLen(i), offsets{i} );
                idl = idl(:,1:newTraceLength);
                [dwt{i},offsets{i}] = idlToDwt(idl);
            end
        end %end for each file.
        
        time = time(1:newTraceLength);

        
    % ---- Extend traces to the same length
    else
        % Extend the traces so they are the same length by duplicating
        % the values in the last frame to pad the end. This can create 
        % problems later on, so use with caution!
        newTraceLength = max( traceLen );

        for i=1:nFiles,
            delta = newTraceLength - numel(data(i).time);
            
            % Resize each data channel in this file.
            for c=1:numel(channelNames),
                chName = channelNames{c};
                data(i).(chName) = [  data(i).(chName) ...
                             repmat( data(i).(chName)(:,end), 1,delta )  ];
            end
        
            % To "extend" the dwell-times, we just have to change the
            % offsets. But this is very easy to do by converting it to an
            % idealization and adding zeros (a marker for regions that are
            % not idealized) to the ends.
            if ~isempty(dwt{i}),
                idl = dwtToIdl( dwt{i}, traceLen(i), offsets{i} );
                idl = [ idl  zeros( size(idl,1), delta )  ];
                [dwt{i},offsets{i}] = idlToDwt(idl);
            end
        end
    end
    
end %if trace length mismatch

waitbar(0.75,h);


% Remove any metadatafields that are not common to all files.
for i=1:nFiles,
    fieldsToRemove = setdiff( fieldnames(data(i).traceMetadata), traceMetadataFields );
    
    if ~isempty(fieldsToRemove),
        data(i).traceMetadata = rmfield( data(i).traceMetadata, fieldsToRemove );
        
        text = sprintf( '%s, ',fieldsToRemove{:} );
        warning( 'combineDatasets:missingMetadataFields', ...
                 'Some metadata fields were removed because they are not in all files: %s', ...
                 text(1:end-2)  );
    end
end

% Merge fluorescence and FRET data and trace metadata. We don't simply copy
% everything from data(1) and overwrite some fields because there may be
% additional fields we don't want (eg, they are not in all fields).
% Warning: that the metadata fields may depend on some fluorescence 
% channels being saved and we don't check for that here!
output.nChannels = numel(channelNames);
output.channelNames = channelNames;
output.time = time; %longest time axis

for c=1:numel(channelNames),
    chName = channelNames{c};
    output.(chName) = vertcat( data.(chName) );
end

if isfield(data,'traceMetadata'),
    output.traceMetadata = [data.traceMetadata];
end

if isfield(data,'fileMetadata'),
    output.fileMetadata = data(1).fileMetadata;
end

clear data;
waitbar(0.9,h);


% Save merged dataset to file.
[p,f,e] = fileparts(outFilename);
if ~isempty(strfind(e,'traces')),
    saveTraces( outFilename, 'traces', output );
elseif ~isempty(strfind(e,'txt')),
    saveTraces( outFilename, 'txt', output );
end


% Merge dwt files if present and consistent.    
n = cellfun( @numel, model ); %count number of states in each model.
if ~all( n(1)==n ),
    warning('.dwt files found for all files, but models have different numbers of states!');

elseif ~all( sampling(1)==sampling ),
    warning('.dwt files found for all files, but they are not the same time resolution and cannot be combined.');

elseif ~all( cellfun(@isempty,dwt) )
    disp('Combining dwt files. I hope you used the same models for these!');

    offsetsAll = [];
    dwtAll     = {};
    modelAll   = zeros( size(model{1}) );

    % Merge all the idealizations, adjusting the offsets.
    fileOffsets = cumsum( [0 nTraces.*newTraceLength] );
    for i=1:numel(offsets),
        offsetsAll = [offsetsAll offsets{i}+fileOffsets(i) ];
        dwtAll     = [dwtAll dwt{i}];
        modelAll   = modelAll + model{i};
    end
    modelAll = modelAll./numel(offsets); %this gives us an "average" model.

    % Save the dwt file.
    if isempty(p), p=pwd; end
    dwtFilename = [p filesep f '.qub.dwt'];
    saveDWT( dwtFilename, dwtAll, offsetsAll, modelAll, sampling(1) );
end


waitbar(1,h);
close(h);



