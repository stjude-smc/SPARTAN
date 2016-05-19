function tdp = tdplot(dwt_input, traces_input, varargin)
% TDPLOT Transition density (TD) plot
%
%   tdp = tdplot(dwtfname,tracefname) creates a transition density plot
%   histogram tdp, which is a two-dimensional contour plot of the FRET value
%   before (x) and after (y) each transition between distinct states.
%   Use tplot(tdp) to display the actual plot.
%
%   tdp = tdplot(IDL, DATA) creates a TD plot using the  state assignment matrix
%   (idealization) IDL and Traces object DATA.
%
%   ... = tdplot(...,options) specifies additional display options as a struct.
% 
%   ... = tdplot(...,'option',value) specifies additional display options as
%   key-value pairs. Defaults are defined in cascadeConstants.
%   
%      normalize determines the units of the output matrix tdp:
%        "total time": transitions per second, counting time before bleaching.
%        "total transitions": percent of all transitions in each histogram bin.
%
%      hideBlinksInTDPlots = true: remove dark state dwells from histogram.
%      truncate_tdplot = true: only consider the first contour_length frames.
%      contour_length: value used above, also shared by makeplots.
%      fret_axis: vector of histogram bin centers.
%      fretField: name of the field in the Traces object to use.
%   
%   See also: tplot, statehist, makeplots.

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% Performance note: textscan > load > dlmread for speed


%% Process input parameters

% Load default values for plotting options.
constants = cascadeConstants;
options = constants.defaultMakeplotsOptions;

% Modify options if they are specified in the argument list.
if nargin==3 && isstruct(varargin{1}),
    options = mergestruct( options, varargin{1} );

elseif nargin>3 && mod(numel(varargin),2)==0,
    options = mergestruct( options, struct(varargin{:}) );
    
elseif nargin>=3
    error('Invalid list of optional parmeters');
end


%% Load FRET and idealizationd data

% Get filenames from user if not passed
if nargin<2    
    dwt_input  = getFile('*.dwt','Choose QuB dwt file:');
    if isempty(dwt_input), return; end
    
    traces_input = getFile('*.traces','Choose an traces file:');
    if isempty(traces_input), return; end
end


% Load FRET data
if ischar(traces_input),
    traces_input = loadTraces(traces_input);
    
elseif ~isa(traces_input,'TracesFret')
    error('tdplot: Invalid traces input');
end

if isfield(options,'fretField') && ~isempty(options.fretField),
    data = traces_input.(options.fretField);
else
    data = traces_input.fret;
end

sampling = traces_input.sampling;
[nTraces,nFrames] = size(data);


% Load dwell-times and convert to a state assignment matrix (idealization).
if ischar(dwt_input),
    [dwt,sampling2,offsets] = loadDWT(dwt_input);
    if sampling~=sampling2,
        warning('Data and idealization time resolution mismatch');
    end
    idl = dwtToIdl(dwt,offsets,nFrames,nTraces);

elseif isnumeric(dwt_input),
    idl = dwt_input;
end


% Optional: truncate the idealization to match the contour plots.
if isfield(options,'truncate_tdplot') && options.truncate_tdplot,
    idl  = idl(:, 1:options.contour_length);
    data = data(:, 1:options.contour_length);
end


%%

% Initialize the 2d histogram
fret_axis = options.fret_axis;

tdp = zeros( numel(fret_axis)+1 );
tdp(1,2:end) = fret_axis;
tdp(2:end,1) = fret_axis;

nTrans = 0;   %total number of transitions
total_time = 0;  %total time in frames


% Count all transitions and add to a transition-density matrix
for i=1:nTraces,
    
    idlTrace = idl(i,:);
    fretData = data(i,:);
    
    % Optional: remove dwells in the dark state, assumed to be state 1.
    if isfield(options,'hideBlinksInTDPlots') && options.hideBlinksInTDPlots,
        fretData = fretData( idlTrace>1 );
        idlTrace = idlTrace( idlTrace>1 );
        if isempty(fretData), continue; end
    end
    
    % Convert back to a dwell-time series.
    dwells = RLEncode( idlTrace );
    states = dwells(:,1);
    times = dwells(states>0,2);
    
    ndwells = numel(times);
    if ndwells==0, continue; end
    
    % Save the number of dwells and total time for later normalization.
    nTrans = nTrans + ndwells-1;
    total_time = total_time + sum(times);    
    
    % Get the start and end indexes of each dwell.
    ti = cumsum( [1 ; times(1:end-1)] );
    tf = cumsum( times );
    
    % For each dwell, calculate the mean FRET value.
    fret = zeros(ndwells,1);
    for j=1:ndwells,
        dwell = fretData(ti(j):tf(j));
        fret(j) = sum(dwell)/numel(dwell);
    end
    
    % Place FRET values into contour bins.
    inds = zeros( ndwells,1 );  %fret bin number assignment of each dwell
    centers = (fret_axis(1:end-1)+fret_axis(2:end)) / 2;

    inds( fret<=centers(1) ) = 1;
    for n=1:numel(centers),
        inds( fret>centers(n) ) = n+1;
    end

    % Add transitions to TD plot
    for k=1:ndwells-1            
        indi = inds(k)+1;    %bin index of FRET data before the transition
        indf = inds(k+1)+1;  %bin index of FRET data after the transition
        tdp(indf,indi) = tdp(indf,indi)+1;
    end

end % for each segment in selection list


% Normalize the plot
tdpData = tdp(2:end,2:end);

if strcmpi(options.normalize,'total time')
    maxVal = total_time*sampling/1000;  %in sec -- independant of framerate

    tdp(1,1) = maxVal;  %save the normalization factor (not plotted)
    tdp(2:end,2:end) = tdpData/double(maxVal);
    
elseif strcmpi(options.normalize,'total transitions')
    maxVal = sum( tdpData(:) )/100; %percent in each bin
    
    tdp(1,1) = maxVal; %max value=total # of transitions
    tdp(2:end,2:end) = tdpData/double(maxVal);

else
    error('Invalid normalization setting');
end



