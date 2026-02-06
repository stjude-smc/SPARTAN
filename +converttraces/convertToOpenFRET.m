function convertToOpenFRET(filenames, varargin)
% convertToOpenFRET  Convert .traces files to OpenFRET format
% This is a wrapper function extracted from TracesToOpenFRET.m to make it
% callable from the App Designer app.
% Author: Zeliha Kilic
% Usage:
%   convertToOpenFRET(filenames, 'title', title, 'description', desc, ...)
% I incorporated out .traces bbinary format within the openfret conversion 
% as described in https://github.com/ajohnsonbuck/openfret-utilities
% See TracesToOpenFRET.m for full documentation of parameters.

convertFiles(filenames, varargin{:});

end

function convertFiles(filenames,varargin) 
% Convert SiMREPS traces.dat files into OpenFRET format
%{
Args 
    Required:
        filenames (cell array of string or char): names of files to be converted to OpenFRET
        path (string or char): path to files

    Optional name-value pairs:
        'compress' (logical): true or false -- specify whether to use zip compression (default = true)
        'title' (str): experiment title
        'description' (str): experiment description
        'experiment_type' (str): type of experiment (e.g., 'SiMREPS')
        'authors' (cell array of str or char): authors of experiment (e.g., {'Jane Doe', 'John Doe'})
        Other recommended metadata:
        'institution' (str): institution where the work was done
        'date' (str): date of experiment
        'experiment_id' (str): Unique ID of experiment
        'buffer_conditions' (str): Buffer used in experiment
        'temperature' (str): Temperature of experiment
        'microscope' (str): Microscope used
        'detector' (str): Detector used
        'objective' (str): Objective used
%}

%% Default attributes
dataset.title = '(Title)';
dataset.description = '(Description of experiment)';
dataset.experiment_type = 'SiMREPS';
channel1.type = 'donor';
channel1.excitation_wavelength = NaN;
channel2.type = 'acceptor';
channel2.excitation_wavelength = NaN;

compress = true; % Compress output to .zip file

%% Parse arguments
if nargin > 2
    for n = 1:2:nargin-2
        if strcmpi(varargin{n},'compress')
            if isa(varargin{n+1},'logical')
               compress = varargin{n+1};
            else
                error("'compress' value must be either true or false")
            end
        elseif strcmpi(varargin{n},'authors') || strcmpi(varargin{n},'author')
            if ~isa(varargin{n+1},'cell')
                dataset.authors = varargin(n+1);
            else
                dataset.authors = varargin{n+1};
            end
        elseif strcmpi(varargin{n},'experiment_id')
            dataset.metadata.experiment_id = varargin{n+1};
        elseif strcmpi(varargin{n},'buffer_conditions')
            dataset.sample_details.buffer_conditions = varargin{n+1};
        elseif strcmpi(varargin{n},'microscope')
            dataset.instrument_details.microscope = varargin{n+1};
        elseif strcmpi(varargin{n},'laser')
            dataset.instrument_details.laser = varargin{n+1};
        elseif strcmpi(varargin{n},'detector')
            dataset.instrument_details.detector = varargin{n+1};
        elseif strcmpi(varargin{n},'objective')
            dataset.instrument_details.objective = varargin{n+1};
        elseif strcmpi(varargin{n},'excitation_wavelength')
            channel1.excitation_wavelength = varargin{n+1}.channel1;
            channel2.excitation_wavelength = varargin{n+1}.channel2;
        elseif strcmpi(varargin{n},'options')
            options = varargin{n+1};
        elseif strcmpi(varargin{n},'title')
            dataset.title = varargin{n+1};
        elseif strcmpi(varargin{n},'description')
            dataset.description = varargin{n+1};
        elseif strcmpi(varargin{n},'date')
            dataset.date = varargin{n+1};
        elseif strcmpi(varargin{n},'institution')
            dataset.institution = varargin{n+1};
        else 
            dataset.(varargin{n}) = varargin{n+1};  % Allow other user-specified fields
        end
    end
end

% Check that at least one channel of data is included; abort operation
% otherwise
if ~(options.use_channel1) && ~(options.use_channel2)
    disp("At least one of options.use_channel1 and options.use_channel2 must be 'true'; aborting operation.");
    return
end


%% Load files and write OpenFRET

traceTemplate = initializeTrace(channel1,channel2);

% Iterate through files, load traces, and write .json for each dataset
for n = 1:numel(filenames)
    fprintf(1,'Creating OpenFRET file for input file %s ...\n',filenames{n});

    filepath = filenames{n};
    if strcmp(filepath(end-3:end),'.dat') || strcmp(filepath(end-3:end),'.mat')
        acceptor = openSiMREPStraces(filepath);
        donor = acceptor*0;
        suffixLength = 3;
    elseif strcmp(filepath(end-6:end),'.traces')
        traces2ch = openTraces(filepath,options.donor_crosstalk);
        donor = traces2ch.donor;
        donorpks = traces2ch.peaks.donor;
        acceptor = traces2ch.acceptor;
        acceptorpks = traces2ch.peaks.acceptor;
        if ~(options.use_channel1)
            donor = donor*0;
            donorpks = donorpks*0;
        elseif ~(options.use_channel2)
            acceptor = acceptor*0;
            acceptorpks = acceptorpks*0;
        end
        suffixLength = 6;
    end
    if all(donor==0)
       channel1.excitation_wavelength = NaN; % Set dummy donor channel wavelength to NaN
    elseif all(acceptor==0)
       channel2.excitation_wavelength = NaN; % Set dummy acceptor channel wavelength to NaN
    end

    ntraces = height(donor);

    dataset.traces = repmat(traceTemplate, 1, ntraces); % Preallocate traces matrix

    fprintf(1,'Loading traces...\n');
    for p = 1:ntraces
        % Place intensity vs time data into channel 1
        dataset.traces(p).channels(1).data = donor(p,:);

        % Place intensity vs time data into channel 2
        dataset.traces(p).channels(2).data = acceptor(p,:);

        % Add peak xy positions, if they are present
        if exist("donorpks","var")
            dataset.traces(p).channels(1).xy = donorpks(p,:);
        end
        if exist("acceptorpks","var")
            dataset.traces(p).channels(2).xy = acceptorpks(p,:);
        end
    end

    % Write to file
    outfilename = strcat(filepath(1:end-suffixLength),'json');
    fprintf(1,'Traces loaded.\nWriting %s to file...\n',outfilename);
    disp(outfilename)
    fid = fopen(outfilename,'w');
    if fid == -1
        error('Cannot write to output directory');
    end
    fclose(fid);
    delete(outfilename);

    write(dataset, outfilename);
    
    % Optionally, compress file
    if compress
        fprintf(1,'Compressing file %s to .zip...\n',outfilename);
        writeZip(outfilename);
    end

    dataset = rmfield(dataset,'traces'); % Reset traces struct for each new file
end

fprintf('Done.\n');

end


function write(dataset, filepath, varargin)
    % Writes a Dataset to a JSON file.
    %
    % Args:
    %   dataset (struct): The Dataset structure.
    %   filepath (str): Path to the JSON file.
    %   optional (str): 'compress' to use .zip compression

    dataset = validateDataset(dataset);
    try
        jsontext = jsonencode(dataset, 'PrettyPrint', true);
        fid = fopen(filepath, 'w');
        fprintf(fid, '%s', jsontext);
        fclose(fid);
        if numel(varargin)>0
            if strcmpi(varargin{1},'compress')
                zipfilename = strcat(filepath,'.zip');
                zip(zipfilename,filepath);
                delete(filepath);
            end
        end
    catch ME
        error('OpenFRET:write:JSONError', 'Error encoding or writing JSON file: %s', ME.message);
    end
end


function dataset = validateDataset(data)
   % Validates Dataset and its components against the schema.
   % This is a rudimentary check, a full schema validation would be ideal.
   if ~isstruct(data) || ~isfield(data, 'title') || ~isfield(data, 'traces')
       error('OpenFRET:validateDataset:InvalidFormat', 'Dataset must be a struct with fields "title" and "traces".');
   end

   for i = 1:length(data.traces)
       trace = data.traces(i);
       if ~isstruct(trace) || ~isfield(trace, 'channels')
            error('OpenFRET:validateDataset:InvalidFormat', 'Trace must be a struct with field "channels".');
       end
       for j = 1:length(trace.channels)
           channel = trace.channels(j);
           if ~isstruct(channel) || ~isfield(channel, 'channel_type') || ~isfield(channel, 'data')
               error('OpenFRET:validateDataset:InvalidFormat', 'Channel must be a struct with fields "channel_type" and "data".');
           end
       end
   end
   dataset = data;
end

function writeZip(filename)
    zipfilename = strcat(filename,'.zip');
    zip(zipfilename,filename);
    delete(filename);
end

function trace = initializeTrace(channel1, channel2)
    % Set properties of trace struct
    trace.channels(1).channel_type = channel1.type;
    trace.channels(1).data = [];
    trace.channels(1).excitation_wavelength = channel1.excitation_wavelength;

    trace.channels(2).channel_type = channel2.type;
    trace.channels(2).data = [];
    trace.channels(2).excitation_wavelength = channel2.excitation_wavelength;
end

function traces = openTraces(filename, donor_crosstalk)
%Args:
    % filename(string or char) = path and filename of .traces file containing two-channel single-molecule data
    %            in binary (ieee-le) format
    % donor_crosstalk(double) = fraction of donor-channel signal that
    %                           appears in acceptor channel (typical value
    %                           = 0.09)

    % Open and read contents from binary traces file

    dataTypes = {'char','uint8','uint16','uint32','uint64', ...
                    'int8', 'int16', 'int32', 'int64', ...
             'single','double','logical','cell','struct'};  %zero-based

    fid = fopen(filename,'r');

    fread(fid,1,'uint32');        % ignore version
    fread(fid,[1,4],'*char');     % ignore signature
    fread(fid,1,'uint16');        % ignore format version

    dataType  = fread( fid, 1, '*uint8' );   % data class (9=single)
    nChannels = fread( fid, 1, '*uint8' );
    nTraces   = fread( fid, 1, 'uint32' );
    traceLen  = fread( fid, 1, 'uint32' );

    szNames      = fread( fid, 1, 'uint32' );
    channelNames = fread( fid, [1,szNames], '*char' );
    channelNames = strsplit(channelNames,char(31), 'CollapseDelimiters',false);

    type = dataTypes{dataType+1};
    data.time = fread( fid, [1,traceLen], type );  %time axis (s)
    
    indexes = 1:nTraces;
    useAll = isequal(indexes, 1:nTraces);

    if dataType==9,
    star = '*';
    else
        disp('there is a problem')
    end

    for i=1:nChannels,
        d = fread( fid, [nTraces,traceLen], [star type] );
        
        if useAll,  %faster if indexing is avoided.
            data.(channelNames{i}) = d;
        else
            data.(channelNames{i}) = d(indexes,:);
        end
    end

    fclose(fid);
    
    % Load peaktable
    try
    filename_pks = strcat(filename(1:end-7), '.pks');
    peaktable = load(filename_pks);
    catch
        fprintf(1, 'Could not load peaktable file named %s; molecule coordinates were not added to dataset.',filename_pks);
    end

    % Store intensities and peak positions in output struct
    traces.donor = data.donor;
    traces.acceptor = data.acceptor;
    if exist("peaktable","var")
        traces.peaks.donor = peaktable(1:2:end,2:3);
        traces.peaks.acceptor = peaktable(2:2:end,2:3);
    else
        traces.peaks.donor = zeros(height(traces.donor),2);
        traces.peaks.acceptor = zeros(height(traces.donor),2);
    end
end

function acceptor = openSiMREPStraces(filepath)
    % Placeholder for SiMREPS traces loading
    % This function would need to be implemented if .dat files are used
    error('SiMREPS .dat file loading not yet implemented in this wrapper');
end

