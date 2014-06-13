function avg_val = avgfret(filename)
%MAKEPLOTS  Creates an array of contour, pop hist, and TD plots
% 
%   AVGFRET(FILENAME)
%   Calculates the average FRET value that would be observed at time=0 in a
%   bulk FRET assay used for drug activity asessment.  FILES specifies the
%   locations of the datasets to plot.  
%   If no FILES are specified, user will be asked to select them.
%   

% TODO: 

%% INITIALIZE & PROCESS FUNCTION PARAMETERS

% If not files specified, prompt user for them.
if ~exist('filename','var'),
    
    [datafile,datapath] = uigetfile({'*.traces'},'Choose a traces file:');
    if datafile==0, return; end  %user hit "cancel"

    filename = fullfile(datapath,datafile);
end

% Strip extensions from files to get the "base name" from which all other
% filenames can be generated (_tdp.txt _hist.txt, etc)
filename = strrep( filename, '.traces', '' );
filename = strrep( filename, '_tdp', '' );
filename = strrep( filename, '.qub', '' );
filename = strrep( filename, '_hist', '' );



%---- USER TUNABLE PARAMETERS ----

% Number of frames summed to give 1D population histogram
pophist_sumlen = 50;

%---------------------------------


% EXTENSIONS for files used:
% dwt_ext  = '.qub.dwt';   % QuB Idealization file
% seg_ext  = '.qub_seg.txt';   % QuB Segments output (whole matrix)
% sel_ext  = '.qub_sel.txt';   % QuB Selection list (not used)
% data_ext = '.qub.txt';   % forQuB raw data (not used)
% raw_ext  = '.txt';   % raw traces data
hist_ext = '_hist.txt';  % TD plot
% tdp_ext  = '.qub_tdp.txt';  % TD plot



%% Make the calculation

% Load population contour plot
hist_filename = [filename hist_ext];

assert( exist(hist_filename,'file')==2, ...
        sprintf('ERROR: Pop histogram file %d does not exist',i) );

cplotdata = load(hist_filename);

%N = sum( cplotdata(2:end,2) );  %number of traces
fretaxis  = cplotdata(2:end,1);
truncated = cplotdata(2:end,2:pophist_sumlen+1);
pophist   = sum(truncated,2);  %/N/pophist_sumlen;

inds = fretaxis>0.05;
pophist = pophist(inds)/sum(pophist(inds));

avg_val = sum( pophist.*fretaxis(inds) );
















