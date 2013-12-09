function constants = cascadeConstants()
% Returns contant used throughput the processing pipeline

constants.version = '2.4';  %pipeline release version number


% ---- Algorithm constants that rarely need to be adjusted.

% Constants for calculating certain properties in traceStat for autotrace.
constants.min_fret = 0.125; % above which consider acceptor dye alive
constants.fretEventTreshold = 0.14; % FRET value used for detecting FRET (binding) events

% acclife (FRET lifetime) does not include short events of FRET. This
% number defines /how/ short to ignore. This is useful to prevent noise
% from contributing too much to acclife. Set to 1 to disable.
constants.rle_min = 5; 

% Number of frames in the window used for calc. background statistics
% in traceStat.m.
constants.NBK = 100; 

% Constants for detecting drops in fluorescence associated with
% photobleaching. NSTD is used for finding big steps associated with single
% fluorophores. The lower threshold in overlap_nstd is used for finding
% multiple photobleaching steps that may be smaller.
% See calcLifetime and traceStat for details.
constants.TAU  = 9;  % Median filter window size
constants.NSTD = 8;  % Step detection threshold in st dev of gradient
constants.overlap_nstd=5;

% Used to calculated a threshold for detecting peaks of fluorescence in
% gettraces.m. In standard deviations of background noise.
constants.gettracesThresholdStd = 8;

constants.blink_nstd=4; % set FRET=0 below threshold (donor is blinking)



% ---- Variable values that are specific to your experimental setup.

% Correction factor for fluorophore brightness and detection efficiency.
%   gamma = apparent acceptor brightness / apparent donor brightness.
% This will depend on your specific setup and the fluorophores used. Not
% necessary to change, but some criteria are less biased if this is set
% correctly.
constants.gamma = 1.0;



% ---- Gettraces default settings

% ADU (arbitrary camera intensity units) to photon conversion factor in
% units of ADU/photon. See camera calibration data sheet. This may depend
% on which digitizer is selected! Check camera documentation.
% If no information is available, comment this line out.
params.photonConversion = 100/3.1;   % 10MHz Evolve 512
%params.photonConversion = 100/2.6;   % 5MHz Evolve 512

% Algorithm settings:
params.don_thresh = 0; %auto
params.overlap_thresh = 2.3;
params.nPixelsToSum   = 4;

% Software alignment settings.
params.alignRotate = 0; %it is slow and not that reliable; disable by default.

params.alignment.theta = -4:0.1:4;
params.alignment.dx    = -5:1:5;
params.alignment.dy    = -5:1:5;
params.alignment.sx    = 1;  %magnification (x)
params.alignment.sy    = 1;  %magnification (y)

% Other options:
params.skipExisting = 0;
params.recursive = 0;
params.quiet = 0;
params.saveLocations = 0;

% Create gettraces parameter profiles for various imaging geometries and
% fluorophores. For the quad-view, the order is UL, UR, LL, LR. These will
% become the items in the drop-down menu at the top of gettraces, with the
% "name" field being the text in the dropdown list.
% See gettraces.m for definitions for these parameters. Acceptable channel
% names include: donor, acceptor, donor2, acceptor2, factor. Factor is a
% fluorescence channel that is not part of any FRET pair.
% FIXME: make it possible to include default alignment settings here,
% probably by giving an alignment filename (.mat).
clear p;

p(1).name        = 'Single-channel (Cy3)';
p(1).geometry    = 1;
p(1).chNames     = {'donor'};
p(1).chDesc      = {'Cy3'};
p(1).wavelengths = 532;
p(1).crosstalk   = 0;
p(1).alignTranslate = 0;

p(2) = p(1);
p(2).name = 'Single-channel (Cy5)';
p(2).wavelengths = 640;

p(3).name        = 'Dual-Cam (Cy3/Cy5)';
p(3).geometry    = 2;
p(3).chNames     = {'donor','acceptor'}; %L/R
p(3).chDesc      = {'Cy3','Cy5'};
p(3).wavelengths = [532 640];
p(3).crosstalk   = 0.12; %donor->acceptor only, 630dcxr
p(3).alignTranslate = 1;  % no problems other than being slow.
% Qinsi's correction for uneven sensitivity of the equipment across the 
% field of view in the acceptor (right) side. Fluorescence intensities are
% at each point are scaled by the amount calculated by the function.
% The function values are listed in the same order as the channels above.
p(3).biasCorrection = {  @(x,y) 1,  ...                        %donor, LHS
                         @(x,y) 0.87854+y*9.45332*10^(-4)  };  %acceptor, RHS

p(4).name        = 'Quad-View (Cy2/Cy3/Cy5/Cy7)';
p(4).geometry    = 3;
p(4).chNames     = {'acceptor','acceptor2','donor','factor'}; %UL/UR/LL/LR
p(4).chDesc      = {'Cy5','Cy7','Cy3','Cy2'};
p(4).wavelengths = [640 730 532 473];
p(4).crosstalk   = 0.26; %donor->acceptor only
p(4).alignTranslate = 0;  %alignment code not fully tested yet!

% Add all of the common settings that do not vary.
fnames = fieldnames(params);
for i=1:numel(fnames),
    [ p.(fnames{i}) ] = deal( params.(fnames{i}) );
end

constants.gettraces_profiles = p;

% Set the default settings profile.
constants.gettraces_defaultProfile = 3;   %Dual-View (Cy3/Cy5)
constants.gettracesDefaultParams = p( constants.gettraces_defaultProfile );



% ---- Default selection criteria for autotrace

criteria.eq_overlap  = 0;       % Remove overlapping molecules
criteria.min_corr    = -1.1;    % 
criteria.max_corr    = 0.5;     % D/A correlation < 0.5
criteria.min_snr     = 8;       % SNR over background
criteria.max_bg      = 70;      % Background noise (in std photons!)
criteria.max_ncross  = 4;       % donor blinking events
criteria.min_acclife = 15;      % FRET lifetime

constants.defaultAutotraceCriteria = criteria;




% ---- Constants for display/plotting functions (makeplots)

% default population FRET contour plot paramters (cplot.m)
options.contour_bin_size = 0.03;     % FRET bin size
options.cplot_scale_factor = 8;      % contour level scaling; depends on bin size
options.contour_length = 50;         % # frames to display in contour plots
options.truncate_statehist = true;   % truncate data used for making statehist/tdplots to
options.truncate_tdplot = false;     %   match the displayed contour plot (contour_length).

% FRET axis range in countor, histogram, and TD plots.
options.fret_axis = -0.1:options.contour_bin_size:1.2;  % bins for histogram calculation.
options.fretRange = [-0.1 1.0];  % range of what is actually displayed.

% Remove X frames from beginning of traces to avoid effect of gain drift.
options.pophist_offset = 0;

% default transition density plot parameters (tdplot.m).
% Adjust tdp_max if spots are saturated (raise it) or not visible (lower it).
options.tdp_max = 0.0025;
options.tdp_fret_axis = options.fret_axis;

% colors for statehist, in order
options.colors = [ 0 0 0 ; 0 0.5 0   ; 1 0 0 ; ...
               0    0.7500    0.7500 ; ...
          0.7500         0    0.7500 ; ...
          0.7500    0.7500         0 ; ...
          0.2500    0.2500    0.2500 ];

% OPTIONS
options.no_statehist = false;  % do not use state occupancy histograms
options.no_tdp       = false;  % do not use TD plots
options.ignoreState0 = true;   % do not include the first (lowest FRET) state in
                               %    state occupancy histograms
options.hideText     = false;  % don't display N=, t/s, etc on plots

options.hideBlinksInTDPlots = false;  % hide transitions to dark state in TD plots

options.saveFiles    = false;  % save histogram txt files for plotting in Origin.

if options.ignoreState0,
   options.colors = options.colors(2:end,:); 
end

constants.defaultMakeplotsOptions = options;



% ---- Other settings
if ispc,
    constants.modelLocation = 'Z:\SharedDocs\Shared QuB\';
else
    constants.modelLocation = '/home/dsterry/data/Daniel/models/';
end

% For MIL (batch kinetics).
warning off MATLAB:maxNumCompThreads:Deprecated
constants.nProcessors = maxNumCompThreads;
