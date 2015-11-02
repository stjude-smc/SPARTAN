function shist=statehist(dwt_input, traces_input, options)
% STATEHIST  state occupancy FRET histogram
%
%   S = STATEHIST( DWT, DATA, options )   where
%   S is a collection (in cols) of FRET histograms for each state, as
%   calculated from the DWT file.  The first col specifies the bins.
%   DWT is the filename of the idealization file from QuB.
%   DATA is the auto.txt filename containing raw Fluorescence/FRET data.
%
%   OPTIONS (optional), can have any of the following fields:
%    - pophist_sumlen:    number of frames to consider when summing
%         histograms. Used to insure population histograms (e.g., 50 frames)
%         match the state histograms by not considering later data.
% 
%    - fret_axis:  specifies FRET histogram bin centers.
%    - fretField:  name of the field in the Traces object to use.

%   Copyright 2007-2015 Cornell University All Rights Reserved.



%% Process input parameters, setting defaults if not provided.
if nargin<3,
    constants = cascadeConstants;
    options = constants.defaultMakeplotsOptions;
end


%% Get filenames from user if not passed
if nargin<2
    dwt_input = getFile('*.dwt','Choose QuB dwt file:');
    if isempty(dwt_input), return; end
    
    traces_input = getFile('*.traces','Choose traces file:');
    if isempty(traces_input), return; end
end


%% Load data
if ischar(traces_input) || isa(traces_input,'TracesFret'),
    % Load Traces object from input or file and select appropriate FRET field.
    if ischar(traces_input),
        traces_input = loadTraces(traces_input);
    end

    if isfield(options,'fretField') && ~isempty(options.fretField),
        fret = traces_input.(options.fretField);
    else
        fret = traces_input.fret;
    end
    
elseif isnumeric(traces_input),
    % FRET data provided directly.
    fret = traces_input;

else
    error('statehist: Invalid traces input');
end

[nTraces,traceLen] = size(fret);


% --- Load the dwell-time data and create an idealization
if ischar(dwt_input),
    [dwt,~,offsets] = loadDWT(dwt_input);
    idl = dwtToIdl(dwt,offsets,traceLen,nTraces);
    
elseif isnumeric(dwt_input),
    idl = dwt_input;

else
    error('statehist: Invalid idealization input');
end
nStates = max(idl(:));


% --- Truncate traces so they match the display length in makeplots
if isfield(options,'truncate_statehist') && options.truncate_statehist,
    traceLen = options.contour_length;
    idl  = idl(:,1:traceLen);
    fret = fret(:,1:traceLen);
end


%% Create normalized histograms for each state.
fret_axis = options.fret_axis;
shist = zeros( numel(fret_axis), nStates+1 );
shist(:,1) = fret_axis;

idl(idl==0&fret~=0) = 1; %otherwise the zero state is just blinks!

for j=1:nStates,
    newdata = hist( fret(idl==j), fret_axis ) /numel(fret);
    shist(:,j+1) = newdata;
end


