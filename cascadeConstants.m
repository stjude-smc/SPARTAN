function constants = cascadeConstants()
% Returns contant used throughput the processing pipeline

constants.version = '2.1.2';  %pipeline release version number


%--- Application settings for memory/CPU usage:

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
constants.photonConversionFactor = 100/3.1; %Fluoescence AU/photon conversion.

constants.overlap_nstd=5; % multiple PB detection threshold
constants.blink_nstd=4;% set FRET=0 below threshold (donor is blinking)


% Donor->Acceptor channel signal crosstalk (CorrectTraces.m)
constants.crosstalk = 0.07;


% CONSTANTS FOR PLOTTING FUNCTIONS

% default transition density plot parameters (tdplot.m)
constants.tdp_fret_axis = -0.1:0.030:1.0;
constants.tdp_max = 0.0025*1;  %note framerate dependance!!

% default population FRET contour plot paramters (cplot.m)
constants.cplot_scale_factor = 8;
constants.contour_length = 50; %default # frames to display in cplot


% 
if ispc,
    constants.modelLocation = 'Z:\SharedDocs\Shared QuB\';
else
    constants.modelLocation = '/home/dsterry/data/Daniel/models/';
end

warning off MATLAB:maxNumCompThreads:Deprecated
constants.nProcessors = maxNumCompThreads;
