function integrateAndSave(this, outname)
%integrateAndSave  Sum single molecule PSFs, save as .traces file.
%
%   integrateAndSave(STKDATA, PARAMS, TRCNAME)
%
%   STKDATA is a fully-initialized MovieParser object, which means that the 
%   openStk, getPeaks, and getIntegrationWindows methods have already been
%   called. The following properties are used for analysis:
%     stkData.movie: Movie object referencing the underlying movie data.
%     stkData.peaks: molecule locations (x first col, y second col).
%     stkData.regionIdx: integration regions, specified as a list of
%                        flat indices for each molecule and each channel.
%     stkData.background: estimated background field image.
%     stkData.bgMask: (optional) specifies background regions to sum to
%                      general a background intensity trace.
%
%   TRCNAME is a the path to output .rawtraces file.
%
%   PARAMS is a structure of optional parameters. See cascadeConstants.m.
%   
%   ALGORITHM:
%   For each location in "peaks", sum a region that includes most of the
%   intensity for that spot and repeat for each time point to get a
%   fluorescence-time trace for each peak. Background subtraction, crosstalk
%   correction, and calculation of derived signals (FRET traces) is all done
%   here. Then the result is saved as a .rawtraces file with metadata.

%   Copyright 2007-2022 All Rights Reserved.


% Process input arguments
narginchk(2,2);
params = this.params;
this.chExtractor.verify();


%% Prepare trace output
quiet = params.quiet;
nTraces = size(this.peaks,1);
nFrames = this.nFrames;  %should be configurable!

% Create channel name list for the final data file. This includes FRET channels,
% which are not in the movie. chNames includes only fluorescence fields.
chNames = this.roles;
nCh = numel(chNames);
dataNames = chNames;  %will include derived traces like fret.

lasers = this.chExtractor.lasers;


%% Integrate fluorescence intensity in movie to create fluorescence traces
if ~quiet,
    wbh = parfor_progressbar(1.1*nFrames,'Starting parallel pool');
end

% Parallelize large movies, where disk access is faster than image processing. 
if nTraces*nFrames/2000 > 1500 && cascadeConstants('enable_parfor')
    pool = gcp;
    M = pool.NumWorkers;
else
    M = 0;  %use GUI thread
end

if ~quiet, wbh.message='Extracting traces from movie data'; end

% The estimated background image is also subtracted to help with molecules
% that do not photobleach during the movie.
% fnames = this.fnames(idxFields);
[bgTrace,bgFieldIdx,bgMask] = deal([]);

if isfield(params,'bgTraceField') && ~isempty(params.bgTraceField),
    bgFieldIdx = params.bgTraceField;
    bgTrace = zeros(nFrames,1,'single');
    bgMask = this.bgMask{bgFieldIdx};
end



% Get a list of field locations (integration windows) to sum into traces.
traces = zeros(nTraces,nFrames,nCh, 'single');
idx = this.regionIdx;  %cell array of channels with [pixel index, molecule id] 

% parfor (k=1:nFrames, M)
for k=1:nFrames
    % Retrieve next frame and separate fluorescence channels.
    frame = this.chExtractor.read(k);
    
    for c=1:nCh
        % Sum intensity within the integration window of each PSF
        traces(:,k,c) = sum( single(frame{c}(idx{c})), 1 );       %#ok<PFBNS>
    end
    
    % Sum intensity from background regions
    if ~isempty(bgFieldIdx)
        bgTrace(k) = mean( single(frame{bgFieldIdx}(bgMask)) );
    end
    
    if mod(k,10)==0 && ~quiet,
        iterate(wbh,10);
    end
end

if ~quiet,
    wbh.message = 'Correcting traces and calculating FRET...';
end

timeAxis = this.chExtractor.movie.timeAxis(1:nFrames);

% Subtract local background
bg = this.chExtractor.background;  %(idxFields);
for c=1:nCh
    bgt = sum( bg{c}(idx{c}), 1);
    traces(:,:,c) = bsxfun(@minus, traces(:,:,c), to_col(bgt) );
end

% Convert from arbitrary camera units to photon counts.
channels = this.chExtractor.channels;
ppc = reshape( [channels.photonsPerCount], [1 1 nCh] );
traces = bsxfun( @times, traces, ppc );

% Handles alternating laser excitation modes
if numel(lasers)>2
    error('Only 2-color ALEX currently supported');
elseif numel(lasers)==2
    % Split trace data by illumination mode
    tracesSplit = cell( numel(lasers), 1 );
    for i=1:numel(lasers)
        idx = lasers(i).framesActive;
        idx = idx(idx<=nFrames);
        tracesSplit{i} = traces(:, idx, i:end);
        
        if i==1, timeAxis=timeAxis(idx); end
    end
    
    % Recreate traces array
    traces = cat(3, tracesSplit{:});
    
    % Update list of channel names for 2-color ALEX.
    nFrames = numel(idx);
    dataNames = [dataNames 'acceptorDirect' 'stoichiometry'];
    chNames   = [chNames   'acceptorDirect'];
end


%% Apply corrections and calculate FRET

% Prompt user for time axis info if not available from movie metadata.
if timeAxis(1)==1 && ~quiet
    disp('Time axis information is not present in movie!');
    a = inputdlg('Time resolution (ms):','No time axis in movie!');
    a = str2double(a);
    if isnan(a)
        disp([mfilename ': Invalid time resolution.']);
    elseif ~isempty(a)
        timeAxis = a*( 0:nFrames-1 );
    end
end

% Create traces object, where the data will be stored.
if ismember('donor',dataNames)
    if ismember('acceptor',dataNames)
        dataNames = [dataNames 'fret'];
    end
    if ismember('acceptor2',dataNames)
        dataNames = [dataNames 'fret2'];
    end
    if any(ismember({'acceptor2','factor'},dataNames))
        data = TracesFret4(nTraces,nFrames,dataNames);
    else
        data = TracesFret(nTraces,nFrames,dataNames);
    end
else
    data = TracesFluor(nTraces,nFrames,dataNames);
end
data.time = timeAxis;

if ~isempty(bgTrace),
    data.fileMetadata.bgTrace = bgTrace;
end

if isfield(params,'zeroMethod'),
    data.fileMetadata.zeroMethod = params.zeroMethod;
end

% Add trace data to the Traces output object
for i=1:numel(chNames),
    data.(chNames{i}) = traces(:,:,i);
end

% Correct for non-uniform camera sensitivity across the field-of-view.
% See the biasCorrection field in cascadeConstants.m for details.
% if isfield(params,'biasCorrection') && ~isempty(params.biasCorrection),
%     for i=1:nCh,
%         chName = data.channelNames{i};
%         x = this.peaks(:,1,i);
%         y = this.peaks(:,2,i);
%         corr = params.biasCorrection{i}(x,y);
%         data.(chName) = data.(chName)  ./  repmat( corr, [1 nFrames] );
%         %data.(chName) = bsxfun( @rdivide, data.(chName), corr );  %TESTME
%     end
% end

% Extract relevant correction parameters from imaging profile
scaleFluor = params.scaleFluor( this.idxActiveChannels );
crosstalk  = params.crosstalk( this.idxActiveChannels, this.idxActiveChannels );

% Subtract background, apply crosstalk/scaling corrections, and calculate FRET.
data = correctTraces( bgsub(data), crosstalk, to_col(scaleFluor) );
data.recalculateFret();



%% Compile trace metadata and save output to disk
if ~quiet,
    wbh.iterate(0.1*nFrames);
    wbh.message = 'Saving traces...';
end

data.fileMetadata.units = 'photons';
data.fileMetadata.wavelengths = [channels.wavelength];
data.fileMetadata.chDesc = {channels.name};
data.fileMetadata.geometry = this.chExtractor.fieldArrangement;  %fixme: ignore=0?
[data.fileMetadata.nY,data.fileMetadata.nX] = size( frame{1} );
data.fileMetadata.profile = params.name;

% Save molecule locations of the picked peaks for later lookup.
for i=1:size(this.peaks,2)
    ch = data.channelNames{i};
    x = num2cell( this.peaks(:,1,i) );
    y = num2cell( this.peaks(:,2,i) );
    [data.traceMetadata.([ch '_x'])] = x{:};
    [data.traceMetadata.([ch '_y'])] = y{:};
end

% Create unique identifiers to track origin of each trace ('movie.tif#123')
for i=1:nTraces
    data.traceMetadata(i).ids = sprintf( '%s#%d', outname, i );
end

% Save data to file.
saveTraces(outname, data);

if ~quiet, close(wbh); end


end %function integrateAndSave


