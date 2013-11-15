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

% Donor->Acceptor channel signal crosstalk (gettraces).
% This applies to both 2-color and 3-color FRET.
constants.crosstalk = 0.07;

% ADU (arbitrary camera intensity units) to photon conversion factor in
% units of ADU/photon. See camera calibration data sheet. This may depend
% on which digitizer is selected! Check camera documentation.
% If no information is available, comment this line out.
constants.photonConversionFactor = 100/3.1;   % 10MHz Evolve 512
%constants.photonConversionFactor = 100/2.6;   % 5MHz Evolve 512



% ---- Gettraces default settings

% Algorithm settings:
params.don_thresh = 0; %auto
params.overlap_thresh = 2.3;
params.nPixelsToSum   = 4;
params.crosstalk = constants.crosstalk;
params.photonConversion = constants.photonConversionFactor;
params.geometry = 2; %dual-channel by default.

% Options for alignment, etc:
params.alignTranslate = 1;  % no problems other than being slow.
params.alignRotate = 0;     % allow field rotation for alignment. not perfect yet.
params.skipExisting = 0;
params.recursive = 0;
params.quiet = 0;
params.saveLocations = 0;

% Channel naming and description. There are specific acceptable names
% (donor, acceptor, etc) that will be used for accessing traces data.
% Descriptions are just metadata and can be anything.
% For now these are just the fluorophores used.
% Dual-channel. Order is: L, R.
constants.gettraces_chNames2 = {'donor','acceptor'};
constants.gettraces_chDesc2  = {'Cy3','Cy5'};
% Three- or Four-channel. Order is: UL, LL, LR, UR.
% Use null ('') for unused channels (three-color).
constants.gettraces_chNames4 = {'donor','factor','acceptor',''};
constants.gettraces_chDesc4 = {'Cy3','Cy2','Cy5',''};

if params.geometry==2,
    params.chNames = constants.gettraces_chNames2;
    params.chDesc  = constants.gettraces_chDesc2;
elseif params.geometry>2,
    params.chNames = constants.gettraces_chNames4;
    params.chDesc  = constants.gettraces_chDesc4;
end

constants.gettracesDefaultParams = params;



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

options.saveFiles    = true;  % save histogram txt files for plotting in Origin.

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
