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
% The first group of settings below will overwrite the ones in the
% profiles, so don't put anything here unless it isn't specified in the
% profile settings!!!

% ADU (arbitrary camera intensity units) to photon conversion factor in
% units of ADU/photon. See camera calibration data sheet. This may depend
% on which digitizer is selected! Check camera documentation.
% If no information is available, comment this line out.
commonParams.photonConversion = 100/3.1;   % 10MHz Evolve 512
%commonParams.photonConversion = 100/2.6;   % 5MHz Evolve 512

% Algorithm settings:
commonParams.don_thresh = 0; %auto
commonParams.overlap_thresh = 2.3;
commonParams.nPixelsToSum   = 4;

commonParams.alignment.theta = -4:0.1:4;
commonParams.alignment.dx    = -4:1:4;
commonParams.alignment.dy    = -4:1:4;  %all dx and dy must be integers
commonParams.alignment.sx    = 1;  %magnification (x)
commonParams.alignment.sy    = 1;  %magnification (y)

% Other options:
commonParams.skipExisting = 0;
commonParams.recursive = 0;
commonParams.quiet = 0;
commonParams.saveLocations = 0;



% Create gettraces parameter profiles for various imaging geometries and
% fluorophores. For the quad-view, the order is UL, UR, LL, LR. These will
% become the items in the drop-down menu at the top of gettraces, with the
% "name" field being the text in the dropdown list.
% See gettraces.m for definitions for these parameters. Acceptable channel
% names include: donor, acceptor, donor2, acceptor2, factor. Factor is a
% fluorescence channel that is not part of any FRET pair.
% chNames has one entry for every subfield (channel) in the field-of-view,
% where some can be empty, meaning to ignore it. The other parameters just
% has values for the channels IN USE.
clear p;
p(1).name        = 'Single-channel (Cy3)';
p(1).geometry    = 1;
p(1).chNames     = {'donor'};
p(1).chDesc      = {'Cy3'};
p(1).wavelengths = 532;
p(1).alignRotate = 0;
p(1).alignTranslate = 0;


p(2) = p(1);
p(2).name        = 'Single-channel (Cy5)';
p(2).wavelengths = 640;
p(2).chDesc      = {'Cy5'};


p(3).name        = 'Dual-Cam (Cy3/Cy5)';
p(3).geometry    = 2;
p(3).chNames     = {'donor','acceptor'}; %L/R
p(3).chDesc      = {'Cy3','Cy5'};
p(3).wavelengths = [532 640];
p(3).crosstalk   = 0.12; %donor->acceptor only
p(3).alignTranslate = 1;  % no problems other than being slow.
p(3).alignRotate = 0;
% Qinsi's correction for uneven sensitivity of the equipment across the 
% field of view in the acceptor (right) side. Fluorescence intensities are
% at each point are scaled by the amount calculated by the function.
% The function values are listed in the same order as the channels above.
p(3).biasCorrection = {  @(x,y) ones(size(x)),  ...            %donor, LHS
                         @(x,y) 0.87854+y*9.45332*10^(-4)  };  %acceptor, RHS

                     
p(4).name        = 'Quad-View (Cy2/Cy3B/Cy5/Cy7)';
p(4).geometry    = 3;
p(4).chNames     = {'acceptor','acceptor2','donor','factor'}; %UL/UR/LL/LR
p(4).chDesc      = {'Cy5','Cy7','Cy3','Cy2'};
p(4).wavelengths = [640 730 532 473];
p(4).crosstalk   = zeros(4);
p(4).crosstalk(3,1) = 0.26;  %Cy3->Cy5
p(4).crosstalk(1,2) = 0.06;   %Cy5->Cy7 (is this correct???)
p(4).alignRotate = 0;
p(4).alignTranslate = 0;


% Add all of the common settings that do not vary.
fnames = fieldnames(commonParams);
for i=1:numel(fnames),
    [ p.(fnames{i}) ] = deal( commonParams.(fnames{i}) );
end


% Set the default settings profile.
constants.gettraces_profiles = p;
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
options.cplot_scale_factor = 7;      % contour level scaling; depends on bin size
options.contour_length = 50;         % # frames to display in contour plots
options.truncate_statehist = true;   % truncate data used for making statehist/tdplots to
options.truncate_tdplot = false;     %   match the displayed contour plot (contour_length).

% This option will normalize the contour plots to the plot with the most
% molecules. This will give a better sense of when there is no FRET in a
% particular condition. The better way to do this is to not have any FRET
% lifetime criteria in autotrace and change the cplot_scale_factor to
% compensate for lower density!
options.cplot_normalize_to_max = false;

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
