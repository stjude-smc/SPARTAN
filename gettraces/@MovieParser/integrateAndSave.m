function integrateAndSave(stkData, outname, params)
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

%   Copyright 2007-2017 Cornell University All Rights Reserved.


% Process input arguments
narginchk(3,3);
if nargin<3
    params = stkData.params;
else
    stkData.params = params;
end



%% Prepare trace output
movie = stkData.movie;
quiet = params.quiet;
nTraces = size(stkData.peaks,1);
nFrames = movie.nFrames;

% Create channel name list for the final data file. This includes FRET channels,
% which are not in the movie. chNames includes only fluorescence fields.
chNames = params.chNames( ~cellfun(@isempty,params.chNames) );
nCh = numel(chNames);
dataNames = chNames;

% Create traces object, where the data will be stored.
if ismember('donor',dataNames)
    if ismember('acceptor',dataNames),
        dataNames = [dataNames 'fret'];
    end
    if ismember('acceptor2',dataNames),
        dataNames = [dataNames 'fret2'];
    end
    
    if params.geometry==4
        data = TracesFret4(nTraces,nFrames,dataNames);
    else
        data = TracesFret(nTraces,nFrames,dataNames);
    end
else
    data = TracesFluor(nTraces,nFrames,dataNames);
end

% If not available from movie, prompt user for time resolution
if ~quiet && movie.timeAxis(1)==1,
    disp('Time axis information is not present in movie!');
    a = inputdlg('Time resolution (ms):','No time axis in movie!');
    a = str2double(a);
    if isnan(a),
        disp('gettraces: Invalid time resolution.');
        return;
    elseif isempty(a), return; %user hit cancel.
    end
    data.time = a*( 0:nFrames-1 );
else
    data.time = movie.timeAxis;
end



%% Integrate fluorescence intensity in movie to create fluorescence traces
if ~quiet,
    wbh = parfor_progressbar(1.1*nFrames,'Extracting traces from movie data');
end

% Parallelize large movies, where disk access is faster than image processing. 
if nTraces*nFrames/2000 > 1500 && cascadeConstants('enable_parfor') && isa(movie,'Movie_TIFF'),
    pool = gcp;
    M = pool.NumWorkers;
else
    M = 0;  %use GUI thread
end

% The estimated background image is also subtracted to help with molecules
% that do not photobleach during the movie.
fnames = stkData.fnames(params.idxFields);
[bgTrace,bgFieldIdx,bgMask] = deal([]);

if isfield(params,'bgTraceField') && ~isempty(params.bgTraceField),
    bgFieldIdx = find( cellfun(@(x)strcmpi(x,params.bgTraceField),fnames), 1,'first' );
    if isempty(bgFieldIdx)
        warning('Background trace field must be an active field');
    else
        bgTrace = zeros(nFrames,1,'single');
        bgMask = stkData.bgMask{bgFieldIdx};
    end
end

traces = zeros(nTraces,nFrames,nCh, stkData.movie.precision);
idx = stkData.regionIdx;  %cell array of channels with [pixel index, molecule id] 

parfor (k=1:nFrames, M)
    for c=1:nCh
        % Extract subfield from this frame and subtract background image
        frame = subfield(movie,fnames{c},k);   %#ok<PFBNS>
        
        % Sum intensity within the integration window of each PSF
        traces(:,k,c) = sum( frame(idx{c}), 1 );       %#ok<PFBNS>
    
        % Sum intensity from background regions
        if ~isempty(bgFieldIdx) && c==bgFieldIdx
            bgTrace(k) = mean( frame(bgMask) );
        end
    end
    
    if mod(k,10)==0 && ~quiet,
        iterate(wbh,10);
    end
end

% Subtract local background
traces = single(traces);
for c=1:nCh
    bgt = sum( stkData.background{c}(idx{c}), 1);
    traces(:,:,c) = bsxfun(@minus, traces(:,:,c), to_col(bgt) );
end


%% Apply corrections and calculate FRET
if ~quiet,
    wbh.message = 'Correcting traces and calculating FRET...';
end

if ~isempty(bgTrace),
    data.fileMetadata.bgTrace = bgTrace;
end

% Convert fluorescence to arbitrary units to photon counts.
if isfield(params,'photonConversion') && ~isempty(params.photonConversion),
    traces = traces./params.photonConversion;
    data.fileMetadata.units = 'photons';
else
    data.fileMetadata.units = 'AU';
end

% Add trace data to the Traces output object
for i=1:nCh,
    data.(chNames{i}) = traces(:,:,i);
end

% Correct for non-uniform camera sensitivity across the field-of-view.
% See the biasCorrection field in cascadeConstants.m for details.
if isfield(params,'biasCorrection') && ~isempty(params.biasCorrection),
    for i=1:nCh,
        chName = data.channelNames{i};
        x = stkData.peaks(:,1,i);
        y = stkData.peaks(:,2,i);
        corr = params.biasCorrection{i}(x,y);
        data.(chName) = data.(chName)  ./  repmat( corr, [1 nFrames] );
        %data.(chName) = bsxfun( @rdivide, data.(chName), corr );  %TESTME
    end
end

if isfield(params,'zeroMethod'),
    data.fileMetadata.zeroMethod = params.zeroMethod;
end

% Subtract background, apply crosstalk/scaling corrections, and calculate FRET.
data = correctTraces( bgsub(data), params.crosstalk, to_col(params.scaleFluor));
data.recalculateFret();



%% Compile trace metadata and save output to disk
if ~quiet,
    wbh.iterate(0.1*nFrames);
    wbh.message = 'Saving traces...';
end

data.fileMetadata.wavelengths = params.wavelengths;
data.fileMetadata.chDesc = params.chDesc;

% Save molecule locations of the picked peaks for later lookup.
for i=1:nCh,
    ch = data.channelNames{i};
    x = num2cell( stkData.peaks(:,1,i) );
    y = num2cell( stkData.peaks(:,2,i) );
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


