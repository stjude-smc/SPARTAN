function shist=statehist(dwtfilename, traces_input, options)
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

%   Copyright 2007-2015 Cornell University All Rights Reserved.



% Get filenames from user if not passed
if nargin<2
    dwtfilename = getFile('*.dwt','Choose QuB dwt file:');
    if isempty(dwtfilename), return; end
    
    tracefilename = getFile('*.traces','Choose traces file:');
    if isempty(tracefilename), return; end
end


% --- Set default parameter values if none provided.
if nargin<3,
    constants = cascadeConstants;
    options = constants.defaultMakeplotsOptions;
end



% --- Open the corresonding qub data file
if isstruct(traces_input)
    data = loadTraces(tracefilename);
    fret = data.fret;  clear data;
    
elseif isnumeric(traces_input),
    fret = traces_input;
    
else
    error('statehist: Invalid traces input');
end

[nTraces,traceLen] = size(fret);



% --- Load the dwell-time data and create an idealization
[dwt,~,offsets] = loadDWT(dwtfilename);
idl = dwtToIdl(dwt,offsets,traceLen,nTraces);
nStates = max(idl(:));


% --- Truncate traces so they match the display length in makeplots
if isfield(options,'contour_length') && isfield(options,'truncate_statehist') && options.truncate_statehist,
    traceLen = options.contour_length;
    idl  = idl(:,1:traceLen);
    fret = fret(:,1:traceLen);
end


% --- Create normalized histograms for each state.
fret_axis = options.fret_axis;
shist = zeros( numel(fret_axis), nStates+1 );
shist(:,1) = fret_axis;

idl(idl==0&fret~=0) = 1; %otherwise the zero state is just blinks!

for j=1:nStates,
    newdata = hist( fret(idl==j), fret_axis ) /numel(fret);
    shist(:,j+1) = newdata;
end

% --- Save resulting histograms.
if nargout==0,
    outfile=strrep(dwtfilename,'.dwt','_shist.txt');
    dlmwrite(outfile,shist,' ');
end


% NOTE: the t here is all transitions, including to 0-FRET
% fprintf('t=%d  N=%d  t/n=%f\n', [t nsegs t/nsegs]);

