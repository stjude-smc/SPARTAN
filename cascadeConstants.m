 function constants = cascadeConstants()
% Returns contant used throughput the processing pipeline

constants.version = '2.2.1';  %pipeline release version number


% ---- Algorithm constants that rarely need to be adjusted.

% Correction factor for fluor detection efficiency.
% 
constants.gamma = 1.0; %

constants.min_fret = 0.125; % above which consider acceptor dye alive
constants.fretEventTreshold = 0.14; % FRET value used for detecting FRET events
constants.rle_min = 5; % count alive when above min_fret for # frames

constants.NBK=100; % Size of window used for calc. background statistics

constants.NSTD=8; % PB detection threshold (see CalcLifetime)
constants.TAU=9; % median filter window size (PB detection).

constants.gettracesThresholdStd = 8; %see gettraces.m

constants.overlap_nstd=5; % multiple PB detection threshold
constants.blink_nstd=4; % set FRET=0 below threshold (donor is blinking)


% ---- Variable values that are specific to your experimental setup.

% Donor->Acceptor channel signal crosstalk (gettraces)
constants.crosstalk = 0.07;

% ADU (arbitrary camera intensity units) to photon conversion. Includes
% dividing by the EM gain (100x). See camera calibration data sheet for ADU
% conversion. This may depend on which digitizer is selected!
% If no information is available, leave this blank 
constants.photonConversionFactor = 100/3.1;   % 10MHz
%constants.photonConversionFactor = 100/2.6;   % 5MHz 



% ---- Gettraces default settings

% Algorithm settings:
params.don_thresh = 0; %auto
params.overlap_thresh = 2.3;
params.nPixelsToSum   = 4;
params.crosstalk = constants.crosstalk;
params.photonConversion = constants.photonConversionFactor;
params.geometry = 2; %dual-channel by default.

% Options for alignment, etc:
params.alignTranslate = 0;
params.alignRotate = 0;
params.refineAlign = 0;
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



% ---- Gettraces default settings
criteria.overlap = 1; % Remove overlapping molecules
criteria.min_corr=-1.1;    % 
criteria.max_corr=0.5;     % D/A correlation < 0.5
criteria.min_snr=8;       % SNR over background
criteria.max_bg=1500;      % Background noise
criteria.max_ncross = 4;   % donor blinking events
criteria.min_acclife = 15; % FRET lifetime

constants.defaultAutotraceCriteria = criteria;



% ---- Constants for display/plotting functions

% default transition density plot parameters (tdplot.m)
constants.tdp_fret_axis = -0.1:0.030:1.0;
constants.tdp_max = 0.0025*1;  %note framerate dependance!!

% default population FRET contour plot paramters (cplot.m)
constants.cplot_scale_factor = 8;
constants.contour_length = 50; %default # frames to display in cplot



% ---- Other settings
if ispc,
    constants.modelLocation = 'Z:\SharedDocs\Shared QuB\';
else
    constants.modelLocation = '/home/dsterry/data/Daniel/models/';
end

warning off MATLAB:maxNumCompThreads:Deprecated
constants.nProcessors = maxNumCompThreads;
