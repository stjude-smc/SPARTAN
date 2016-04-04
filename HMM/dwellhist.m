function varargout = dwellhist(dwtfilename, inputParams)
%dwellhist  Dwell-time histograms
% 
%   [X,HIST] = dwellhist(FILES) dwell-time histograms for each .dwt file in 
%   the cell array FILES (must all have same time resolution). HIST is a cell
%   array with files in rows and states in columns size=[nFiles,nStates].
%   The X are the bin edges in seconds.
%
%   dwellhist() prompts the user for a list of files.
%
%   dwellhist(...,PARAMS) give optional parameters in struct PARAMS:
%
%      'removeBlinks': Remove dwells in dark states (class 1). Default=true.
%                      Dwells broken up by such blinks are merged, with the
%                      time during the blink added to surrounding dwells.
%
%      'logX':         Use a log-scale times axis to aid visualization of 
%                      multi- exponential distributions (default=true).
%                      See Sigworth and Sine (1987), Biophys J 50, p. 1047-1054.
%
%      'dx':           Log-scale time axis bin size. default=0.2.
%
%      'normalize':    log-scale histogram normalization method:
%                      'off'   - raw dwell counts, no normalization.
%                      'state' - each state; sum of each histogram=1. (default)
%                      'file'  - all states; sum of all histograms per file=1.
%                      'time'  - dwell counts per second of observation time
%
%   See also: dwellplots, lifetime_exp, loadDwelltimes, removeBlinks.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


%% ---- USER TUNABLE PARAMETERS ----

params.logX = true;
params.dx = 0.2;  %log-scale bin width (0.1=25%, 0.2=60%, 0.5=3-fold, 1=10-fold)
params.removeBlinks = true;
params.normalize = 'state';

% Merge options, giving the user's options precedence.
if nargin>1,
    params = mergestruct( params, inputParams );
end

% Check parameters
assert( all(ismember(params.normalize,{'off','state','file','time'})), ...
        'Invalid normalization option' );



%% Prompt user for file names if not given.

% Get list of .dwt files to load
if nargin<1,
    dwtfilename = getFiles('*.dwt','Choose dwell-time files');
end
if ischar(dwtfilename), dwtfilename={dwtfilename}; end
dwtfilename = findDwt(dwtfilename);
if numel(dwtfilename)==0, return; end


% If there are no outputs requested, display instead.
if nargout==0,
    dwellplots(dwtfilename,params);
end



%%
nFiles = numel(dwtfilename);
dwells  = cell(nFiles,1);  %consolidated list of dwell times in each state
sampling = zeros(nFiles,1);

for i=1:nFiles,
    % Load dwell-times a list per state, concatinating dwells from all traces,
    % ignoring the zero-FRET state if applicable.
    if params.removeBlinks,
        [dwells{i},sampling(i)] = loadDwelltimes( dwtfilename{i}, 'removeBlinks' );
        dwells{i} = dwells{i}(2:end);
    else
        [dwells{i},sampling(i)] = loadDwelltimes( dwtfilename{i} );
    end
end

% Verify all input idealizations are roughly consistent.
if ~all(sampling==sampling(1)),
    warning('Mismatched time resolution of input files');
end
sampling = sampling(1)/1000;  %convert to seconds.

nStates = cellfun(@numel, dwells);
if ~all(nStates==nStates(1)),
    warning('Idealization models have a different number of states');
end
nStates = max(nStates);


% Get dwell time limits for setting axes limits later.
maxTime = 0;  %longest dwell in seconds
totalTime = zeros(nFiles,nStates);

for i=1:nFiles,
    dwellc = dwells{i};
    maxTime = max( maxTime, max(vertcat(dwellc{:})) );
    totalTime(i,:) = cellfun(@sum, dwellc)';
end



%% Calculate dwell time bins (EDGES)
if ~params.logX,
    % Linear X-axis in seconds.
    dwellaxis = 0:sampling:maxTime;
else
    % Create a log time axis with a fixed number of bins.
    % histcounts uses bin edges of E(k) <= X(i) < E(k+1).
    dwellaxis = log10(sampling):params.dx:log10(maxTime*3);
    
    % Force the bins edges to be exact intervals of the time resolution.
    % The histogram will better sample the discrete nature of the data.
    maxFrames = ceil(maxTime*3/sampling);
    fullaxis = log10( (1:maxFrames)*sampling )';
    dwellaxis = unique( nearestBin(dwellaxis, fullaxis) );
%     dwellaxis = unique( floor(fullaxis/sampling)*sampling );
    
    % Normalization factor to account for varying-sized bins.
    dlx = dwellaxis(2:end) - dwellaxis(1:end-1);
    dlx = [dlx dlx(end)];
end


%% Calculate histograms
histograms = cell(nFiles,nStates);

for file=1:nFiles,
    ndwells = cellfun(@numel,dwells{file});
            
    for state=1:nStates,
        % Small constant ensures dwells fall in the correct histogram bin.
        dwellc = dwells{file}{state} +sampling/10;
        
        % Make linear-scale survival plot.
        if ~params.logX,
            counts = histc( dwellc, dwellaxis );
            histdata = sum(counts) - cumsum(counts);
            histdata = histdata/histdata(1);
        
        % Make log-scale Sine-Sigworth plot (linear ordinate).
        else
            counts = histc( log10(dwellc)', dwellaxis );
            histdata = counts./dlx;  %normalize by log-space bin size
            histdata = histdata/sum(histdata);  %normalize to 1
            
            switch params.normalize
                case 'off'  %raw dwell counts
                    histdata = histdata*ndwells(state);
                case 'state'  %fraction of counts in each bin for this state
                    histdata = 100*histdata;
                case 'file'  %fraction of counts in each bin across entire file
                    histdata = 100*histdata *ndwells(state)/sum(ndwells);
                case 'time'  %fraction of dwells per total observation time
                    histdata = histdata*ndwells(state)/sum(totalTime(file,:));
            end
        end
        
        histograms{file,state} = to_col(histdata);
    end
end


% Combine histograms into a matrix for saving.
if params.logX,
    dwellaxis = 10.^dwellaxis;
end

output = {dwellaxis,histograms};
[varargout{1:nargout}] = output{1:nargout};


end %function dwellhist



function [newVal,idx] = nearestBin( values, bins )
% For each VALUE, find the BIN with the closest value.

newVal = zeros( size(values) );
idx = zeros( size(values) );

for i=1:numel(values),
    [~,idx(i)] = min( abs(bins-values(i)) );
    newVal(i) = bins(idx(i));
end

end



