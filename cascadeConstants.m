function constants = cascadeConstants()
% Returns contant used throughput the processing pipeline

constants.version = '2.7';  %pipeline release version number



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



% ---- Variable values that are specific to your experimental setup.

% Correction factor for fluorophore brightness and detection efficiency.
%   gamma = apparent acceptor brightness / apparent donor brightness.
% This will depend on your specific setup and the fluorophores used. Not
% necessary to change, but some criteria are less biased if this is set
% correctly.
constants.gamma = 1.0;





%% ======================  Gettraces Default Settings ====================== %%
% The first group of settings below will overwrite the ones in the
% profiles, so don't put anything here unless it isn't specified in the
% profile settings!!!

% Placeholders
commonParams = struct( 'name','', 'geometry',0, 'idxFields',[], 'chNames',{}, ...
                       'chDesc',{}, 'wavelengths',[], 'crosstalk',[], ...
                       'biasCorrection',{} );

% ADU (arbitrary camera intensity units) to photon conversion factor in
% units of ADU/photon. See camera calibration data sheet. This may depend
% on which digitizer is selected! Check camera documentation.
% If no information is available, comment this line out.
commonParams(1).photonConversion = 100/3.1;   % 10MHz Evolve 512
%commonParams(1).photonConversion = 100/2.6;   % 5MHz Evolve 512

% Algorithm settings:
commonParams.don_thresh     = 0;   %molecule detection threshold (0=automatic)
commonParams.overlap_thresh = 2.3; %remove molecules that are w/i X pixels.
commonParams.nPixelsToSum   = 4;   %number of pixels to sum per trace
commonParams.alignMethod    = 1;   %disabled, assume aligned.

% Other options:
commonParams.skipExisting  = 0; %batch mode: skip files already processed.
commonParams.recursive     = 0; %batch mode: search recursively.
commonParams.quiet         = 0; %don't output debug messages.
commonParams.saveLocations = 0; %save molecule locations to a text file


% Create gettraces parameter profiles for various imaging geometries and
% fluorophores. For the quad-view, the order is UL, UR, LL, LR. These will
% become the items in the drop-down menu at the top of gettraces, with the
% "name" field being the text in the dropdown list.
% See gettraces.m for definitions for these parameters. Acceptable channel
% names include: donor, acceptor, donor2, acceptor2, factor. Factor is a
% fluorescence channel that is not part of any FRET pair.
% idxFields specifies the mapping between channel names and the physical
% position on the CCD chip. For the Quad-View the order is UL/UR/LL/LR.
% For all parametes, only list channels that will be used.
% NOTE: you must put the channels in their spectral order.
clear p; clear profiles


p = commonParams;
p.name        = 'EMCCD, Single-channel (Cy3)';
p.geometry    = 1;
p.idxFields   = 1; %only one channel
p.chNames     = {'donor'};
p.chDesc      = {'Cy3'};
p.wavelengths = 532;
profiles(1) = p;


p.name        = 'EMCCD, Single-channel (Cy5)';
p.wavelengths = 640;
p.chDesc      = {'Cy5'};
profiles(end+1) = p;


p = commonParams;
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


p.name        = 'EMCCD 5MHz, Dual-Cam (Cy3/Cy5)';
p.photonConversion = 100/2.6;   % 5MHz Evolve 512, assuming 100x gain.
profiles(end+1) = p;


p = commonParams;
p.name        = 'Quad-View (Cy2/Cy3/Cy5/Cy7)';
p.geometry    = 3;
p.idxFields   = [4 3 1 2];  %field order: LR/LL/UL/UR
p.chNames     = {'donor','acceptor','donor2','acceptor2'};
p.chDesc      = {'Cy2','Cy3','Cy5','Cy7'};
p.wavelengths = [473 532 640 730];
p.crosstalk   = zeros(4);
p.crosstalk(2,3) = 0.13;   %Cy3->Cy5
p.crosstalk(3,4) = 0.06;   %Cy5->Cy7 (is this correct???)
profiles(end+1) = p;


p = commonParams;
p.name        = 'Quad-View (Cy3/Cy5 only)';
p.geometry    = 3;
p.idxFields   = [3 1]; %field order: LL/UL
p.chNames     = {'donor','acceptor'};
p.chDesc      = {'Cy3','Cy5'};
p.wavelengths = [532 640];
p.crosstalk   = 0.13;   %Cy3->Cy5
profiles(end+1) = p;


p = commonParams;
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


p = commonParams;
p.name        = 'Quad-View (Cy5/Cy7)';
p.geometry    = 3;
p.idxFields   = [1 2]; % field order: LL/UL/LL
p.chNames     = {'donor','acceptor'};
p.chDesc      = {'Cy5','Cy7'};
p.wavelengths = [640 730];
p.crosstalk   = 0.11;
profiles(end+1) = p;


p = commonParams;
p.name        = 'sCMOS, Twin-Cam (Cy3/Cy5)';
p.geometry    = 2;
p.idxFields   = [1 2]; %L/R
p.chNames     = {'donor','acceptor'};
p.chDesc      = {'Cy3','Cy5'};
p.wavelengths = [532 640];
p.crosstalk   = 0.115;  %donor->acceptor
p.nPixelsToSum = 5;  %5-6
p.photonConversion = 2.04; %0.49 e-/ADU
profiles(end+1) = p;

% This last one is for a few movies we took in the beginning, where the cameras
% were reversed in metamorph (Cy5 on the left). Not needed anymore?
p.name        = 'sCMOS, Twin-Cam (Cy3/Cy5) REVERSED';
p.idxFields   = [2 1]; %L/R; they're reversed. oops.
profiles(end+1) = p;


% Set the default settings profile.
constants.gettraces_profiles = profiles;
constants.gettraces_defaultProfile = 3;   %Dual-View (Cy3/Cy5)
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
warning off MATLAB:maxNumCompThreads:Deprecated
constants.nProcessors = maxNumCompThreads;






