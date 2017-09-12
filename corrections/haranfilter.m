function data = haranfilter( input, M,P,windowSizes )
% haranfilter  Filter noise in fluorescence traces, preserving anticorrelation.
%
%  DATA = haranfilter( DATA )
%  applies a non-linear filter specifically designed to reduce noise in
%  single-molecule fluorescence traces, while preserving sharp transitions
%  in FRET due to conformational transitions. Uses default parameter values.
%  FIXME: this synthax no longer works!
%
%  DATA = haranfilter( FILES )
%  filters the fluorescence data for each file in FILES cell array. The output
%  is saved as a new file with the extension "_flt.traces".
%
%  DATA = haranfilter( ..., M, P, win )
%  K and M are window sizes (larger values average over larger regions,
%  good for high-noise data with few transitions). P is an amplification
%  factor -- the higher the value, the more transitions will be accentuated
%  as sharp, anticorrelated changes in fluorescence.
%  Defaults: M=10, P=5.
%
% See G. Haran, Chem Phys 307 (204): 137-145. "Noise reduction in
% single-molecule fluorescence trajectories of folding proteins".

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Setup default parameter values.
% FIXME: default values from the paper, but not sure what would work best.
if nargin<2,
    M = 10;
end
if nargin<3,
    P = 5;
end
if nargin<4,
    windowSizes = [4 8 16 32 64]; %window sizes to use
end

% Verify input parameters
assert( numel(M)==1&numel(P)==1, 'Parameters M and P must be scalars' );
assert( M>1&P>0, 'M and P must be positive numbers' );


% If no files given, obtain a list from the user.
if nargin<1,
    input = getFiles;
    if numel(input)==0, return; end
end

if nargin>=1 && isa(input,'Traces'),
    % Handle trace data being given instead of a set of files.
    data = copy(input);
    wbh = parfor_progressbar( data.nTraces, 'Filtering fluorescence data...' );
    data = haranfilter_file( data, M,P,windowSizes, wbh );
    close(wbh)
    return;
    
elseif ischar(input),
    files = {input};
    
elseif iscell(input),
    files = input;
else
    error('Invalid input');
end


% Get the total number of traces so we can calibrate the waitbar.
nTraces = sizeTraces(files);
wbh = parfor_progressbar( sum(nTraces), 'Filtering fluorescence data...' );

try
    % For each file, filter fluorescence intensity & save result.
    for i=1:numel(files),
        % Load fluorescence data
        data = loadTraces( files{i} );

        % Filter data
        data = haranfilter_file(data, M,P,windowSizes, wbh);

        % Save resulting data.
        [p,f,ext] = fileparts( files{i} );
        outFilename = fullfile( p, [f '_flt' ext]);
        saveTraces( outFilename, data );
    end
catch e
    if strcmpi(e.identifier,'parfor_progressbar:cancelled')
        disp( [mfilename ': Operation cancelled by user'] );
        data = [];
    else
        rethrow(e);
    end
end

close(wbh);

end





%%
function data = haranfilter_file( data, M,P,windowSizes, wbh )
% Wrapper function to run the filter on a single file.

% Start parallel pool for large data sets, if enabled.
constants = cascadeConstants;
if data.nTraces*data.nFrames/2000 > 5 && constants.enable_parfor,
    pool = gcp;
    nProc = pool.NumWorkers;
else
    nProc = 0;
end

% Filter the data.
d = data.donor;
a = data.acceptor;

parfor (i=1:data.nTraces, nProc)
% for i=1:data.nTraces,
    [d(i,:),a(i,:)] = haranfilter3( d(i,:),a(i,:), M,P,windowSizes );
    wbh.iterate();  %update waitbar. consider using mod().
end
data.donor = d;
data.acceptor = a;

% Calculate FRET from the filtered fluorescence traces.
% FIXME: The noise now has strange properties, so recalculateFret doesn't work.
% data.recalculateFret();
fret = data.fret;
data.fret = a./(a+d);
data.fret( fret==0 ) = 0;

end



%%
function [Id_out,Ia_out] = haranfilter3( Id,Ia, M,P,windowSizes )
% This is the actual Haran filter that removes noise in donor and acceptor
% traces between anticorrelated changes.
% FIXME: this can probably be optimized!
    
% Pad the data with zeros so filtering can be performed on the edges.
Id = [zeros(size(Id,1),M) Id];
Ia = [zeros(size(Ia,1),M) Ia];

[~,traceLen] = size(Id);
K = numel(windowSizes);

%% Calculate forward and backward predictors (weighted averages).
Ifd = zeros( K,traceLen ); %forward  predictors (donor dye)
Ibd = zeros( K,traceLen ); %backward predictors (donor dye)
Ifa = zeros( K,traceLen ); %forward  predictors (acceptor dye)
Iba = zeros( K,traceLen ); %backward predictors (acceptor dye)

for k=1:K,  %for each window size...
    N = windowSizes(k);
    
    for i=1+N:traceLen-N,  %for each point in the original data...
        Ifd(k,i) = sum( Id(i-N:i-1) )/N;
        Ibd(k,i) = sum( Id(i+1:i+N) )/N;
        
        Ifa(k,i) = sum( Ia(i-N:i-1) )/N;
        Iba(k,i) = sum( Ia(i+1:i+N) )/N;
    end
    
end


%% Calculate forward and backward weights.
wf = zeros( K,traceLen ); %forward predictor weights
wb = zeros( K,traceLen ); %backward predictor weights

for k=1:K,  %for each window size...    
    for i=1+M:traceLen-M,
        j = 0:M-1;
        wf(k,i) = sum(  ( Id(i-j)-Ifd(k,i-j) ).^2 + ( Ia(i-j)-Ifa(k,i-j) ).^2 ); %slow!
        wb(k,i) = sum(  ( Id(i+j)-Ibd(k,i+j) ).^2 + ( Ia(i+j)-Iba(k,i+j) ).^2 ); %slow!
    end
end

wf = wf.^(-P);
wb = wb.^(-P);


% Normalize predictor weights over all window sizes.
% sum weights over k (rows), giving one normalization factor for each datapoint (columns).
C = sum( wf+wb ); 
C = repmat( C, [K,1] );

wf = wf./C;
wb = wb./C;


%% Calculate filtered donor and acceptor traces
% Weighted average of forward and backward predictors.
Id_out = sum( wf.*Ifd + wb.*Ibd );
Ia_out = sum( wf.*Ifa + wb.*Iba );


% Remove padding from edges...
Id_out = Id_out(:,M+1:end);
Ia_out = Ia_out(:,M+1:end);


end

