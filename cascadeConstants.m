function constants = cascadeConstants()
% Returns contant used throughput the processing pipeline

constants.version = '2.9';  %pipeline release version number


% FIXME: this script is called many times and is getting large enough to add
% time to some function calls. Consider using a static version of the return
% parameters. Problem: result is not updated when the file is changed.


%% ======================  Global Algorithm Settings ====================== %%
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




%% ======================  Gettraces Default Settings ====================== %%
% Create gettraces parameter profiles for various imaging geometries and
% fluorophores. Some parameters, especially crosstalk, may depend on the
% specific biological system, fluorophores, and filters used. Others, like
% photonConversion, nPixelsToSum, nHoodSize, and overlap_thresh depend on the
% cameras used, level of binning, magnification, etc.
% 
% Channel names must be listed in spectral order (blue, green, red, IR).
% The order these channels appear in the movie is specified by the "idxFields"
% property. The order is: UL, UR, LL, LR.   (For two fields: L, R).
% For example a value of [2 4 1] specifies gives a field order of: UR,LR,UL.
% 
% See gettraces.m for definitions for these parameters. Acceptable channel
% names include: donor, acceptor, acceptor2, and factor.
% This allows for 2- or 3-color FRET (but not four-color with two FRET pairs).
% Factor is a fluorescence channel that is not part of any FRET pair.
%

%------------------------
% Default settings are given here. Unless another value is given in the profile
% definition, these values are used.

cmosCommon = struct( 'name','', 'geometry',0, 'idxFields',[], 'chNames',{}, ...
                       'chDesc',{}, 'wavelengths',[], 'crosstalk',[], ...
                       'biasCorrection',{} );

% Gettraces GUI settings:
cmosCommon(1).alignMethod = 1;  %disabled, assume aligned.
cmosCommon.skipExisting   = 0;  %batch mode: skip files already processed.
cmosCommon.recursive      = 0;  %batch mode: search recursively.
cmosCommon.quiet          = 0;  %don't output debug messages.
cmosCommon.saveLocations  = 0;  %save molecule locations to a text file

% Conversion from camera units (ADU) to photons (photoelectrons).
% See camera calibration datasheet. May depend on which digitizer is selected!
% If no information is available, comment this line out.
cmosCommon.photonConversion = 2.04;   % 0.49 e-/ADU  (manual says 0.46?)

% Algorithm settings:
% These depend on the PSF size relative to pixel size and must be optimized.
cmosCommon.don_thresh     = 0;   %molecule detection threshold (0=automatic)
cmosCommon.overlap_thresh = 3.5; %remove molecules that are w/i X pixels.
cmosCommon.nPixelsToSum   = 9;   %number of pixels to sum per trace
cmosCommon.nhoodSize      = 2;   %integrate within this neighborhood (px distance from peak)
                                   %  1=3x3 area, 2=5x5 area, 3=7x7 area, etc.


% Default settings for EMCCD (Evolve 512) cameras with 2x2 binning.
emccdCommon = cmosCommon;
emccdCommon.photonConversion = 100/3.1;  % 10MHz chipset, gain 4 (3x), 100x EM gain
emccdCommon.overlap_thresh   = 2.3;      % 
emccdCommon.nPixelsToSum     = 4;        % optimal SNR
emccdCommon.nhoodSize        = 1;        % 3x3 area



%------------------------
% Settings for particular setups are listed here. Each entry is concatinated
% onto the end of the list (profiles).

clear p; clear profiles

%------  sCMOS cameras  -------
p = cmosCommon;
p.name         = 'sCMOS, Single-channel (Cy3)';
p.geometry     = 1;
p.idxFields    = 1; %only one channel
p.chNames      = {'donor'};
p.chDesc       = {'Cy3'};
p.wavelengths  = 532;
profiles(1)    = p;


p.name        = 'sCMOS, Twin-Cam (Cy3/Cy5)';
p.geometry    = 2;
p.idxFields   = [1 2]; %L/R
p.chNames     = {'donor','acceptor'};
p.chDesc      = {'Cy3','Cy5'};
p.wavelengths = [532 640];
p.crosstalk   = 0.115;  %donor->acceptor (no bandpass filters!)
profiles(end+1) = p;

% % For a few movies taken with old versions of Flash Gordon
% p.name        = 'sCMOS, Twin-Cam (Cy3/Cy5) REVERSED';
% p.idxFields   = [2 1]; %R/L
% profiles(end+1) = p;


p.name        = 'sCMOS, Multi-Cam (Cy3/Cy5/Cy7, NO bandpass)';
p.geometry    = 3;
p.idxFields   = [1 2 4]; % field order: UL,UR,LR.
p.chNames     = {'donor','acceptor','acceptor2'};
p.chDesc      = {'Cy3','Cy5','Cy7'};
p.wavelengths = [532 640 730];
p.crosstalk   = zeros(4);
p.crosstalk(1,2) = 0.11;   %Cy3->Cy5 (same as 2-color; was 0.066 with bandpasses in?)
p.crosstalk(2,3) = 0.04;   %Cy5->Cy7 (0.015 with bandpasses in?)
profiles(end+1) = p;

% % For a few movies taken with old versions of Flash Gordon
% p.name        = 'sCMOS, Multi-Cam (Cy3/Cy5/Cy7) OLD ORDER';
% p.idxFields   = [3 2 1]; % field order: LL, UR, UL
% profiles(end+1) = p;



%------  EMCCD cameras  ------
p = emccdCommon;
p.name        = 'EMCCD, Single-channel (Cy3)';
p.geometry    = 1;
p.idxFields   = 1; %only one channel
p.chNames     = {'donor'};
p.chDesc      = {'Cy3'};
p.wavelengths = 532;
profiles(end+1) = p;


p.name        = 'EMCCD, Single-channel (Cy5)';
p.wavelengths = 640;
p.chDesc      = {'Cy5'};
profiles(end+1) = p;


p = emccdCommon;
p.name        = 'EMCCD, Dual-Cam (Cy3/Cy5)';
p.geometry    = 2;
p.idxFields   = [1 2]; %L/R
p.chNames     = {'donor','acceptor'};
p.chDesc      = {'Cy3','Cy5'};
p.wavelengths = [532 640];
p.crosstalk   = 0.075;  %donor->acceptor
% Qinsi's correction for uneven sensitivity of the equipment across the 
% field of view in the acceptor (right) side. Fluorescence intensities are
% at each point are scaled by the amount calculated by the function.
% The function values are listed in the same order as the channels above.
p.biasCorrection = {  @(x,y) ones(size(x)),  ...            %donor, LHS
                      @(x,y) 0.87854+y*9.45332*10^(-4)  };  %acceptor, RHS
profiles(end+1) = p;


p.name        = 'EMCCD, Dual-Cam (Cy3/Cy5, no binning)';
p.nhoodSize   = 2; %5x5 area
p.nPixelsToSum = 7;
profiles(end+1) = p;


p = emccdCommon;
p.name        = 'Quad-View (Cy3/Cy5 only)';
p.geometry    = 3;
p.idxFields   = [3 1]; %field order: LL/UL
p.chNames     = {'donor','acceptor'};
p.chDesc      = {'Cy3','Cy5'};
p.wavelengths = [532 640];
p.crosstalk   = 0.13;   %Cy3->Cy5
profiles(end+1) = p;


p = emccdCommon;
p.name        = 'Quad-View (Cy3/Cy5/Cy7)';
p.geometry    = 3;
p.idxFields   = [3 1 2]; % field order: LL/UL/LL
p.chNames     = {'donor','acceptor','acceptor2'};
p.chDesc      = {'Cy3','Cy5','Cy7'};
p.wavelengths = [532 640 730];
p.crosstalk   = zeros(4);
p.crosstalk(1,2) = 0.12;   %Cy3->Cy5
p.crosstalk(2,3) = 0.06;   %Cy5->Cy7 (is this correct???)
profiles(end+1) = p;


p.name = 'Quad-View (Cy3/Cy5/Cy7, no binning)';
p.nhoodSize = 2; %5x5 area
p.nPixelsToSum = 7;
profiles(end+1) = p;


p = emccdCommon;
p.name        = 'Quad-View (Cy5/Cy7)';
p.geometry    = 3;
p.idxFields   = [1 2]; % field order: LL/UL/LL
p.chNames     = {'donor','acceptor'};
p.chDesc      = {'Cy5','Cy7'};
p.wavelengths = [640 730];
p.crosstalk   = 0.11;
profiles(end+1) = p;


% Set the default settings profile.
constants.gettraces_profiles = profiles;
constants.gettraces_defaultProfile = 2;   %sCMOS Cy3/Cy5
constants.gettracesDefaultParams = profiles( constants.gettraces_defaultProfile );





%% =================  Autotrace Default Selection Criteria ================= %%

criteria.eq_overlap  = 0;       % Remove overlapping molecules
criteria.min_corr    = -1.1;    % 
criteria.max_corr    = 0.5;     % D/A correlation < 0.5
criteria.min_snr     = 8;       % SNR over background
criteria.max_bg      = 70;      % Background noise (in std photons!)
criteria.max_ncross  = 4;       % donor blinking events
criteria.min_acclife = 15;      % FRET lifetime

constants.defaultAutotraceCriteria = criteria;





%% ======================  Makeplots Default Settings ====================== %%
% ---- Constants for display/plotting functions (makeplots)

% default population FRET contour plot paramters (cplot.m)
options.contour_bin_size = 0.03;     % FRET bin size
options.cplot_scale_factor = 8;      % contour level scaling; depends on bin size
options.contour_length = 50;         % # frames to display in contour plots
options.truncate_statehist = true;   % truncate data used for making statehist/tdplots to
options.truncate_tdplot = false;     %   match the displayed contour plot (contour_length).

% This option will normalize all histograms to the dataset with the most
% molecules. This will give a better sense of when there is no FRET in a
% particular condition. The better way to do this is to not have any FRET
% lifetime criteria in autotrace and change the cplot_scale_factor to
% compensate for lower density!
options.cplot_normalize_to_max = false;  % false by default.

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

if options.ignoreState0,
   options.colors = options.colors(2:end,:); 
end

constants.defaultMakeplotsOptions = options;





%% ============================  Other Settings ============================ %%

if ispc,
    constants.modelLocation = 'Z:\SharedDocs\Shared QuB\';
else
    constants.modelLocation = '/media/Z/SharedDocs/Shared QuB/';
end

% For MIL (batch kinetics).
constants.nProcessors = feature('numCores');

% Set to false to disable parfor, which is slow on some older computers.
constants.enable_parfor = constants.nProcessors>1;






