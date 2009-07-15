function constants = cascadeConstants()
% Returns contant used throughput the processing pipeline

constants.version = '1.5';  %pipeline release version number


% Correction for fluor detection efficiency.
% 
constants.gamma = 0.7; %

constants.min_fret = 0.125; % above which consider acceptor dye alive
constants.fretEventTreshold = 0.14; % FRET value used for detecting FRET events
constants.rle_min = 5; % count alive when above min_fret for # frames

constants.NBK=100; % Size of window used for calc. background statistics

constants.NSTD=8; % PB detection threshold (see CalcLifetime)
constants.TAU=9; % median filter window size (PB detection).

constants.overlap_nstd=5; % multiple PB detection threshold $
constants.blink_nstd=4;% set FRET=0 below threshold (donor is blinking)


% Donor->Acceptor channel signal crosstalk (CorrectTraces.m)
constants.crosstalk = 0.075;


% CONSTANTS FOR PLOTTING FUNCTIONS

% default transition density plot parameters (tdplot.m)
constants.tdp_fret_axis = -0.1:0.030:1.0;
constants.tdp_max = 0.0025*1;  %note framerate dependance!!

% default population FRET contour plot paramters (cplot.m)
constants.cplot_scale_factor = 8;
constants.contour_length = 50; %default # frames to display in cplot



constants.modelLocation = '/home/dsterry/cornell/data/Daniel/models/';
constants.binaryLocation = '/home/dsterry/cornell/code/cascade_binary/';
