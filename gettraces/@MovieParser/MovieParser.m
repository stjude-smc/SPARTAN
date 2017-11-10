classdef MovieParser < handle
% MovieParser   Display an interactive list of fret traces.
%
%   p = MovieParser(FILE) opens the given fluorescence movie file (.tif or .stk)
%   and prepares a parser object that will assist in movie analsis tasks, such
%   as summing fluorescence fields, detecting molecule locations, finding
%   integration windows for each, and saving the traces to file.
%
%   See also: gettraces, gettraces_gui, etc.

%   Copyright 2017 Cornell University All Rights Reserved.



properties (GetAccess=public, SetAccess=protected)
    stage = 0;           % 1=openStk, 2=getPeaks, 3=getIntegrationWindows
    
    % Basic data available when movie is first loaded in openStk()
    nChannels;           % Number of fluorescence channels??
    movie;               % Movie_TIFF or Movie_STK object
    stk_top;             % Sum of the first 10 frames (fields separate).
    background;          % Estimated background image from first 10 frames
    endBackground;       % Last 10 frames, total fluorescence intensity
    
    % Picked molecules from getPeaks()
    total_t;             % Registered, total intensity image used for picking
    peaks;               % Molecule locations in all fields
    total_peaks;         % Molecule locations in the total intensity image
    rejectedPicks;       %
    rejectedTotalPicks;  % Locations of molecules with overlapping PSFs
    fractionOverlapped;  % Fraction rejected of total
    alignStatus;         % Alignment strct: dx, dy, theta, sx, sy, abs_dev, tform, quality
    
    % Integration windows from getIntegrationWindows()
    regionIdx;              % X,Y coordinates of each integration window
    integrationEfficiency;  % Estimated fraction of intensity collected
    fractionWinOverlap;     % Fraction of pixels used by multiple molecules
    bgMask;                 % Logical mask of pixels used for summing background fluorescence
end

properties (GetAccess=public, SetAccess=public)
    params;              % Analysis settings. See cascadeConstants.m
end



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


