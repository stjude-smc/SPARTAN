function shist=statehist(dwtfilename, tracefilename, options)
% STATEHIST  state occupancy FRET histogram
%
%   S = STATEHIST( DWT, DATA, options )   where
%   S is a collection (in cols) of FRET histograms for each state, as
%   calculated from the DWT file.  The first col specifies the bins.
%   DWT is the filename of the idealization file from QuB.
%   DATA is the auto.txt filename containing raw Fluorescence/FRET data.
%   
%
%   OPTIONS (optional), can have any of the following fields:
%    - pophist_sumlen:    number of frames to consider when summing
%         histograms. Used to insure population histograms (e.g., 50 frames)
%         match the state histograms by not considering later data.
% 
%    - fret_axis:  specifies FRET histogram bin centers.
%

% Created: Mar 6, 2008    Daniel Terry



% Get filenames from user if not passed
if ~exist('tracefilename','var')
     %---Open the QuB dwt file from idealization
    [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file:');
    if dwtfile==0, return;  end
        
    dwtfilename=strcat(dwtpath,dwtfile);

    %---Open the corresonding qub data file
    [tracefile tracepath]=uigetfile('*.txt','Choose qub data file:');
    if tracefile==0, return;  end
    
    tracefilename=strcat(tracepath,tracefile);
end


% --- Histogram axes
if nargin<4,
    options = struct();
end

% FIXME: this should be in cascadeConstants.
if  ~isfield(options,'fret_axis'),
    BIN=0.030;
    fret_axis = -0.2:BIN:1.2;
else
    fret_axis = options.fret_axis;
end


% --- Open the corresonding qub data file
data = loadTraces(tracefilename);
fret = data.fret;  clear data;
[nTraces,traceLen] = size(fret);


% --- Load the dwell-time data and create an idealization
[dwt,sampling,offsets,model] = loadDWT(dwtfilename);
idl = dwtToIdl(dwt,traceLen,offsets);

if iscell(model),  model = model{1};  end
nStates = numel(model)/2;


% --- Truncate traces so they match the display length in makeplots
if isfield(options,'pophist_sumlen'),
    traceLen = options.pophist_sumlen;
    idl  = idl(:,1:traceLen);
    fret = fret(:,1:traceLen);
end


% --- Create normalized histograms for each state.
shist = zeros( numel(fret_axis), nStates+1 );
shist(:,1) = fret_axis;

for j=1:nStates,
    newdata = hist( fret(idl==j), fret_axis ) /numel(fret);
    shist(:,j+1) = reshape(newdata, numel(fret_axis), 1);
end

% --- Save resulting histograms.
outfile=strrep(dwtfilename,'.dwt','_shist.txt');
dlmwrite(outfile,shist,' ');


% NOTE: the t here is all transitions, including to 0-FRET
% fprintf('t=%d  N=%d  t/n=%f\n', [t nsegs t/nsegs]);

