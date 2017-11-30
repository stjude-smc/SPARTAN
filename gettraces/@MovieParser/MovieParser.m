classdef MovieParser < handle
% MovieParser   Display an interactive list of fret traces.
%
%   p = MovieParser(FILE) opens the given fluorescence movie file (.tif or .stk)
%   and prepares a parser object that will assist in movie analsis tasks, such
%   as summing fluorescence fields, detecting molecule locations, finding
%   integration windows for each, and saving the traces to file.
%
%   See also: gettraces_gui.

%   Copyright 2017 Cornell University All Rights Reserved.



properties (GetAccess=public, SetAccess=protected)
    % Basic data available when movie is first loaded in openStk().
    % stk_top, backround, stdbg include ALL fields, even unused ones,
    % and are arranged in column-major order (not wavelength order!)
    movie;               % Movie_TIFF or Movie_STK object
    stk_top;             % Sum of the first 10 frames (cell array of fields).
    background;          % Estimated background image from first 10 frames  (cell array of fields)
    stdbg;               % stdev of background noise at the end of movie (last 10 frames)
    nFrames;             % number of actual time units
    
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

% properties (GetAccess=public, SetAccess=public)
%     params;              % Analysis settings. See cascadeConstants.m
% end



methods
    % Constructor
    function this = MovieParser(input, params)
        this = openStk(this,input, params);
    end
    
    % Open movie and prepare for viewing
    this = openStk(this, input, params);
    
    % Detect peaks of intensity in registered, total intensity image
    this = getPeaks(this, params);
    
    % Find integration windows for each molecule
    this = getIntegrationWindows(this, params);
    
    % Sum fluorescence in integration windows and save fluorescence traces
    integrateAndSave(this, filename, params);
end



end %class MovieParser


