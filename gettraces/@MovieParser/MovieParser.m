classdef MovieParser < handle
% MovieParser   Display an interactive list of fret traces.
%
%   p = MovieParser(FILE) opens the given fluorescence movie file (.tif or .stk)
%   and prepares a parser object that will assist in movie analsis tasks, such
%   as summing fluorescence fields, detecting molecule locations, finding
%   integration windows for each, and saving the traces to file.
%
%   See also: Movie, ChannelExtractor, gettraces_gui.

%   Copyright 2017-2022 All Rights Reserved.



properties (GetAccess=public, SetAccess=protected)
    % Basic data available when movie is first loaded in openStk().
    chExtractor;         % ChannelExtractor encapsulates Movie and splits frame data into channels
    stk_top;             % Average of first 10 frames (cell array of fields).
    background;          % Estimated background image from first 10 frames  (cell array of fields)
    stdbg;               % stdev of background noise at the end of movie
    
    % Picked molecules from getPeaks()
    total_t;             % Registered, total intensity image used for picking
    peaks;               % Molecule locations in all fields (molID,dim,channel)
    total_peaks;         % ... in total intensity image
    rejectedPicks;       % Locations of molecules with overlapping PSFs
    rejectedTotalPicks;  % ... in total intensity image
    fractionOverlapped;  % Fraction of molecules rejected due to overlapping PSFs
    alignStatus;         % Alignment struct: dx, dy, theta, sx, sy, abs_dev, tform, quality
    
    % Integration windows from getIntegrationWindows()
    regionIdx;           % X,Y coordinates of each integration window
    bgMask;              % Logical mask of pixels used for summing background fluorescence
    
    % Integration startistics for display, also from getIntegrationWindows()
    integrationEfficiency;  % Estimated fraction of intensity collected
    psfWidth;               % Average number of pixels to integrate 70% of total intensity
    fractionWinOverlap;     % Fraction of pixels used by multiple molecules
end


% These can be directly manipulated by the user (gettraces_gui.m)
% FIXME: use set method to update structure and raise event if changed.
properties (GetAccess=public, SetAccess=public)
    % Correction parameters
    crosstalk;   % Spectral crosstalk fraction (NxN matrix)
    scaling;     % Scale each channel to account of uneven brightness/detection effiency.
    %ade;        % Acceptor direct excitation corrections
    
    %See chExtractor.channels for name,wavelength,photonsPerCount,description,role.
    %These are not carried over across movies (resets each time!)
    
    % General analysis settings. See cascadeConstants.m.
    % These are intended to be carried between movies as long as they're
    % recorded with the same microscope.
    params;
end


properties (Dependent)
    % Number of frames in the movie after deinterlacing colors.
    % This will be equal to movie.nFrames unless fluorescence channels
    % appear as interlaced frames instead of field areas (not common).
    nFrames;
    
    % Index of each channel into params.
    % The values in params include all known cameras for a microscope,
    % not just the ones in the current file.
    idxParams;
end



methods
    % Constructor
    function this = MovieParser(varargin)
        if nargin>0
            this = this.openStk( varargin{:} );
        end
    end
    
    % Open movie and prepare for viewing
    this = openStk(this, input, params);
    
    % Detect peaks of intensity in registered, total intensity image
    this = getPeaks(this);
    
    % Find integration windows for each molecule
    this = getIntegrationWindows(this);
    
    % Sum fluorescence in integration windows and save fluorescence traces
    integrateAndSave(this, filename);
    
    % get/set methods
    function value = get.nFrames(this)
        % Accounts for channels stacked as separate frames
        value = this.chExtractor.nFrames;
    end
    
    function value = get.idxParams(this)
        value = cellfun( @(x)find(strcmpi(x,this.params.chNames)), {this.channels.name} );
        %fixme: use closest wavelengths instead as a fallback.
    end
    
    
    %TODO: Need some new methods to allow user to change geometry and
    %update channel values...
end



end %class MovieParser


