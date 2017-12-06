function output = cascadeConstants(fname)
%cascadeConstants  Parameter values used throughout SPARTAN.
%
%   const = cascadeConstants() is a struct of parameters.

%   Copyright 2007-2017 Cornell University All Rights Reserved.


% Cached for faster execution. Automatically reset if file is modified.
persistent constants;
if isempty(constants),
    constants = makeConstants();
end

% If requested, grab just a specific field
if nargin>0
    output = constants.(fname);
else
    output = constants;
end

end %function cascadeConstants




function constants = makeConstants()

% Random number will change whenever the file is updated
constants.tstamp = now();

% Version info displayed in title bars
constants.version = '3.6';
constants.software = ['Cornell SPARTAN ' constants.version];



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
constants.blink_nstd=4; % set FRET=0 below threshold (donor is blinking)




%% ======================  Gettraces Default Settings ====================== %%
% Definitions of gettraces parameter profiles.
%
% Channels must be listed in spectral order (UV, blue, green, red, IR).
% "geometry" defines how the fluorescence channels are arranged, encoded in the
% size of the matrix, and their wavelength order. For example, [2,3; 1 0] says
% the fields are arranged as four quadrants in each frame, with the following
% wavelength order: lower-left, upper-left, upper-right (lower-right ignored).
% 
% Valid channel names are defined in TracesFret4.m and include:
% donor, acceptor, acceptor2, and factor.

%------------------------
% Default settings are given here. Unless another value is given in the profile
% definition, these values are used.
cmosCommon = struct( 'name','', 'geometry',0, 'idxFields',[], 'chNames',{}, ...
                       'chDesc',{}, 'wavelengths',[], 'crosstalk',[], ...
                       'scaleFluor',[],'biasCorrection',{} );

% Gettraces GUI settings:
cmosCommon(1).alignMethod = 3;  %auto/ICP.
cmosCommon.skipExisting   = 1;  %batch mode: skip files already processed.
cmosCommon.recursive      = 1;  %batch mode: search recursively.
cmosCommon.quiet          = 0;  %don't output debug messages.
cmosCommon.zeroMethod     = 'threshold';  %method for detecting donor blinks
cmosCommon.bgTraceField   = ''; %get a background intensity trace for this field (L,R,TR,etc)

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
cmosCommon.bgBlurSize     = 6;   %background estimation window size (see openStk.m)
cmosCommon.thresh_std     = 8;   %stdev's of background noise above which to pick molecules
cmosCommon.alignment = struct([]);

% Default settings for EMCCD (Evolve 512) cameras with 2x2 binning.
% emccdCommon = cmosCommon;
% emccdCommon.photonConversion = 100/3.1;  % 10MHz chipset, gain 4 (3x), 100x EM gain
% emccdCommon.overlap_thresh   = 2.3;      % 
% emccdCommon.nPixelsToSum     = 4;        % optimal SNR
% emccdCommon.nhoodSize        = 1;        % 3x3 area
% emccdCommon.bgBlurSize       = 6;



%------------------------
% Settings for particular setups are listed here. Each entry is concatinated
% onto the end of the list (profiles).
clear p; clear profiles

%------  sCMOS cameras  -------
p = cmosCommon;
p.name         = 'sCMOS, Single-channel (Cy3)';
p.geometry     = 1;
p.chNames      = {'donor'};
p.chDesc       = {'Cy3'};
p.wavelengths  = 532;
p.scaleFluor   = 1;
profiles(1)    = p;


p.name        = 'sCMOS, Twin-Cam (Cy3/Cy5, L/R)';
p.geometry    = [1 2];
p.chNames     = {'donor','acceptor'};
p.chDesc      = {'Cy3','Cy5'};
p.wavelengths = [532 640];
p.crosstalk   = zeros(2);
p.crosstalk(1,2) = 0.115;  %donor->acceptor (no bandpass filters!)
p.scaleFluor  = [1 1];
profiles(end+1) = p;

p.name        = 'sCMOS, Twin-Cam (Cy3/Cy5, sequential)';
p.geometry    = cat(3,1,2);
profiles(end+1) = p;


p.name        = 'sCMOS, Multi-Cam (Cy2/3/5)';
p.geometry    = [2 3; 1 0];  % field order: LL,UL,UR.
p.chNames     = {'factor','donor','acceptor'};
p.chDesc      = {'Cy2','Cy3','Cy5'};
p.wavelengths = [473 532 640];
p.crosstalk   = zeros(3);
% p.crosstalk(1,2) = 0.0;   %Cy2->Cy3 (FIXME)
p.crosstalk(2,3) = 0.11;   %Cy3->Cy5 (FIXME)
p.scaleFluor  = [1 1 1];
profiles(end+1) = p;


p.name        = 'sCMOS, Multi-Cam (Cy3/Cy5/Cy7)';
p.geometry    = [1 2; 0 3];  % field order: UL,UR,LL.
p.chNames     = {'donor','acceptor','acceptor2'};
p.chDesc      = {'Cy3','Cy5','Cy7'};
p.wavelengths = [532 640 730];
p.crosstalk   = zeros(3);
p.crosstalk(1,2) = 0.11;   %Cy3->Cy5 (same as 2-color; was 0.066 with bandpasses in?)
p.crosstalk(2,3) = 0.04;   %Cy5->Cy7 (0.015 with bandpasses in?)
p.scaleFluor  = [1 1 7];  %Multiply Cy5 and Cy7 traces by this amount.
profiles(end+1) = p;



% %------  EMCCD cameras  ------
% p = emccdCommon;
% p.name        = 'EMCCD, Single-channel (Cy3)';
% p.geometry    = 1;
% p.idxFields   = 1; %only one channel
% p.chNames     = {'donor'};
% p.chDesc      = {'Cy3'};
% p.wavelengths = 532;
% p.scaleFluor  = 1;
% profiles(end+1) = p;
% 
% 
% p.name        = 'EMCCD, Single-channel (Cy5)';
% p.wavelengths = 640;
% p.chDesc      = {'Cy5'};
% profiles(end+1) = p;
% 
% 
% p = emccdCommon;
% p.name        = 'EMCCD, Dual-Cam (Cy3/Cy5)';
% p.geometry    = 2;
% p.idxFields   = [1 2]; %L/R
% p.chNames     = {'donor','acceptor'};
% p.chDesc      = {'Cy3','Cy5'};
% p.wavelengths = [532 640];
% p.crosstalk   = zeros(2);
% p.crosstalk(1,2) = 0.075;  %donor->acceptor
% p.scaleFluor  = [1 1];
% % Qinsi's correction for uneven sensitivity of the equipment across the 
% % field of view in the acceptor (right) side. Fluorescence intensities are
% % at each point are scaled by the amount calculated by the function.
% % The function values are listed in the same order as the channels above.
% p.biasCorrection = {  @(x,y) ones(size(x)),  ...            %donor, LHS
%                       @(x,y) 0.87854+y*9.45332*10^(-4)  };  %acceptor, RHS
% profiles(end+1) = p;
% 
% 
% p.name        = 'EMCCD, Dual-Cam (Cy3/Cy5, no binning)';
% p.nhoodSize   = 2; %5x5 area
% p.nPixelsToSum = 7;
% profiles(end+1) = p;
% 
% 
% p = emccdCommon;
% p.name        = 'Quad-View (Cy3/Cy5 only)';
% p.geometry    = 4;
% p.idxFields   = [3 1]; %field order: LL/UL
% p.chNames     = {'donor','acceptor'};
% p.chDesc      = {'Cy3','Cy5'};
% p.wavelengths = [532 640];
% p.crosstalk   = zeros(2);
% p.crosstalk(1,2) = 0.13;   %Cy3->Cy5
% p.scaleFluor  = [1 1];
% profiles(end+1) = p;
% 
% 
% p = emccdCommon;
% p.name        = 'Quad-View (Cy3/Cy5/Cy7)';
% p.geometry    = 4;
% p.idxFields   = [3 1 2]; % field order: LL/UL/LL
% p.chNames     = {'donor','acceptor','acceptor2'};
% p.chDesc      = {'Cy3','Cy5','Cy7'};
% p.wavelengths = [532 640 730];
% p.crosstalk   = zeros(4);
% p.crosstalk(1,2) = 0.12;   %Cy3->Cy5
% p.crosstalk(2,3) = 0.06;   %Cy5->Cy7 (is this correct???)
% p.scaleFluor  = [1 1 7];
% profiles(end+1) = p;
% 
% 
% p.name = 'Quad-View (Cy3/Cy5/Cy7, no binning)';
% p.nhoodSize = 2; %5x5 area
% p.nPixelsToSum = 7;
% profiles(end+1) = p;
% 
% 
% p = emccdCommon;
% p.name        = 'Quad-View (Cy5/Cy7)';
% p.geometry    = 4;
% p.idxFields   = [1 2]; % field order: LL/UL/LL
% p.chNames     = {'donor','acceptor'};
% p.chDesc      = {'Cy5','Cy7'};
% p.wavelengths = [640 730];
% p.crosstalk   = zeros(2);
% p.crosstalk(1,2) = 0.11;
% p.scaleFluor  = [1 1];
% profiles(end+1) = p;


% Save indexes for getting wavelength-sorted list of selected fields.
% DO NOT EDIT
for i=1:numel(profiles)
    [val,idx] = sort( profiles(i).geometry(:) );
    profiles(i).idxFields = idx(val>0);
end


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

% Name of the field to use for FRET data. (could instead be fret2, etc)
options.fretField = 'fret';

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

% This option will remove zero-FRET states (acceptor or donor dark states) from
% contour plots. This makes it easier to see any changes in state occupancy over
% time without interference from loss of high-FRET over time.
% NOTE that toward the end of the histogram, very few traces may contribute.
options.cplot_remove_bleached = false;

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

% FRET contour plot colormap (cplot and tplot)
options.cmap = [
251 251 196
251 251 196
253 253 232
255 255 255
186 217 254
 74 138 255
  0   0 255
  0 255   0
 71 254  16
255 255   0
255 180  40
255 111  15
255   0   0
255   0   0 ]/255;

% OPTIONS
options.no_statehist = false;  % do not use state occupancy histograms
options.no_tdp       = false;  % do not use TD plots
options.ignoreState0 = true;   % do not include the first (lowest FRET) state in
                               %    state occupancy histograms
options.hideText     = false;  % don't display N=, t/s, etc on plots

options.hideBlinksInTDPlots = false;  % hide transitions to dark state in TD plots
options.normalize = 'total time';  % TD plot units are transitions per second

if options.ignoreState0,
   options.colors = options.colors(2:end,:); 
end

constants.defaultMakeplotsOptions = options;





%% ============================  Other Settings ============================ %%

if isdeployed,
    constants.modelLocation = pwd;
else
    % Blanchard lab model locations.
    % FIXME: this should be removed in the future.
    if ispc,
        constants.modelLocation = 'Z:\SharedDocs\Shared QuB\';
    else
        constants.modelLocation = '/media/Z/SharedDocs/Shared QuB/';
    end
end

% For MIL (batch kinetics).
constants.nProcessors = feature('numCores');

% Set to false to disable parfor, which is slow on some older computers.
constants.enable_parfor = constants.nProcessors>1 & ~isdeployed;


end  %function makeConstants


