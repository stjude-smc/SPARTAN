function gettraces(varargin)
% GETTRACES  Extract smFluorescence traces from movies

%   Copyright 2007-2017 Cornell University All Rights Reserved.

% If calling directly from command line, launch the GUI.
if nargin==0,
    gettraces_gui;
else
    error('Interface is depricated. Use MovieParser');
end

