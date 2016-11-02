function integrateAndSave(stkData, stk_fname, params)
%integrateAndSave  Sum single molecule PSFs, save as .traces file.
%
%   integrateAndSave(STK, PEAKS, FILENAME, PARAMS)
%
% For each location in "peaks", sum a region that includes most of the
% intensity for that spot and repeat for each time point to get a
% fluorescence-time trace for each peak. Background subtraction, crosstalk
% correction, and calculation of derived signals (FRET traces) is all done
% here. Then the result is saved as a .rawtraces file with metadata.

%   Copyright 2007-2016 Cornell University All Rights Reserved.

constants = cascadeConstants;
movie = stkData.movie;
nFrames = movie.nFrames;

% Start the progress bar before initial setup; indicate something is happening.
quiet = params.quiet;
if ~quiet,
    wbh = parfor_progressbar(1.1*nFrames,'Extracting traces from movie data');
end

% Get x,y coordinates of picked peaks
peaks = stkData.peaks;
Npeaks = size(peaks,1);
x = peaks(:,1);
y = peaks(:,2);


% Create channel name list for the final data file. This includes FRET channels,
% which are not in the movie. chNames includes only fluorescence fields.
chNames = params.chNames( ~cellfun(@isempty,params.chNames) );
nCh = numel(chNames);
nTraces = Npeaks/nCh;

dataNames = chNames;

if ismember('acceptor',dataNames),
    dataNames = [dataNames 'fret'];
end

if ismember('acceptor2',dataNames),
    dataNames = [dataNames 'fret2'];
end

% Create traces object, where the data will be stored.
if params.geometry==4 && ismember('donor',dataNames)
    data = TracesFret4(nTraces,nFrames,dataNames);
elseif ismember('donor',dataNames)
    data = TracesFret(nTraces,nFrames,dataNames);
else
    data = TracesFluor(nTraces,nFrames,dataNames);
end

data.time = movie.timeAxis;

if isfield(params,'zeroMethod'),
    data.fileMetadata(1).zeroMethod = params.zeroMethod;
end

% Ask the user
if ~quiet && data.time(1)==1,
    disp('Time axis information is not present in movie!');
    a = inputdlg('Time resolution (ms):','No time axis in movie!');
    a = str2double(a);
    if ~isempty(a), %empty if user hit Cancel or entered an invalid number
        data.time = a*( 0:nFrames-1 );
    else
        disp('Using frames as the time axis');
    end
end

% Parallelize very large TIFF movies, where disk access is quick compared to
% image processing. 
if nTraces*nFrames/2000 > 1500 && constants.enable_parfor && isa(movie,'Movie_TIFF'),
    % Processing large TIFF movies is CPU limited. Use parfor to parallelize.
    pool = gcp;
    M = pool.NumWorkers;
else
    % For small datasets, do everything in the GUI thread (regular for loop).
    M = 0;
end

% Create a trace for each molecule across the entire movie.
% The estimated background image is also subtracted to help with molecules
% that do not photobleach during the movie.
traces = zeros(Npeaks,nFrames,'single');

doBgTrace = isfield(params,'bgTraceField') && ~isempty(params.bgTraceField);
if doBgTrace,
    bgTrace = zeros(nFrames,1,'single');
    bgMask = imerode(stkData.bgMask, ones(3));  %avoid PSF tails
    bgMask = bgMask & subfield_mask(bgMask, params.bgTraceField);
else
    bgMask = [];
end

idx = stkData.regionIdx;  %pixel, peak(chId:nCh:end).
bg = single(stkData.background);
nPx = params.nPixelsToSum;

parfor (k=1:nFrames, M)
    % NOTE: 25% faster by converting to int16, with no change to sCMOS data.
    % But EMCCD have slight differences due to 15-bit overflows?
    frame = single(movie.readFrame(k)) - bg; %#ok<PFBNS>
    
    if nPx>1,
        traces(:,k) = sum( frame(idx) );
    else
        traces(:,k) = frame(idx);
    end
    
    if doBgTrace,
        bgTrace(k) = mean( frame(bgMask) );
    end
    
    % Update waitbar. Using mod speeds up the loop, but isn't ideal because
    % indexes are executed somewhat randomly. Reasonably accurate despite this.
    if mod(k,10)==0 && ~quiet,
        wbh.iterate(10); %#ok<PFBNS>
    end
end
if ~quiet,
    wbh.message = 'Correcting traces and calculating FRET...';
end

if doBgTrace,
    data.fileMetadata.bgTrace = bgTrace;
end

% Convert fluorescence to arbitrary units to photon counts.
if isfield(params,'photonConversion') && ~isempty(params.photonConversion),
    traces = traces./params.photonConversion;
    data.fileMetadata(1).units = 'photons';
else
    data.fileMetadata(1).units = 'AU';
end


% Add trace data to the Traces output object
for i=1:nCh,
    data.(chNames{i}) = traces(i:nCh:end,:);
end


% Correct for non-uniform sensitivity across the field-of-view in the donor
% and acceptor fields. This is generally fixed and determined by the
% optical properties in the light path and sometimes the CCD chips. The
% parameters in cascadeConstants (.biasCorrection) are lamda functions that
% give a scaling factor for each point in the field-of-view. To generate
% this, find functions that give a flat response profile. The elements in
% the lamba function cell array are in the same order as the channels.
if isfield(params,'biasCorrection') && ~isempty(params.biasCorrection),
    for i=1:nCh,
        chName = data.channelNames{i};
        corr = params.biasCorrection{i}( x(i:nCh:end), y(i:nCh:end) );
        data.(chName) = data.(chName)  ./  repmat( corr, [1 nFrames] );
    end
end

% Subtract background, apply crosstalk/scaling corrections, and calculate FRET.
data = bgsub(data);
data = correctTraces(data, params.crosstalk, to_col(params.scaleFluor));
data.recalculateFret();


% ---- Metadata: save various metadata parameters from movie here.
if ~quiet,
    wbh.iterate(0.1*nFrames);
    wbh.message = 'Saving traces...';
end

data.fileMetadata.wavelengths = params.wavelengths;
data.fileMetadata.chDesc = params.chDesc;


% -- Save the locations of the picked peaks for later lookup.
tempx = num2cell(x);
tempy = num2cell(y);

for i=1:nCh,
    ch = data.channelNames{i};
    [data.traceMetadata.([ch '_x'])] = tempx{i:nCh:end};
    [data.traceMetadata.([ch '_y'])] = tempy{i:nCh:end};
end


% -- Create unique trace identifiers.
% This is just the full path to the movie plus a trace number. This can be
% used to later find the corresponding original movie data for each
% individual trace, even after many rounds of processing.
for i=1:size(data.donor,1),
    data.traceMetadata(i).ids = sprintf( '%s#%d', stk_fname, i );
end

% ---- Save data to file.
[p,name]=fileparts(stk_fname);
save_fname = fullfile(p, [name '.rawtraces']);
saveTraces(save_fname, data);

if ~quiet,
    close(wbh);
end
% disp(toc);

end %function integrateAndSave
