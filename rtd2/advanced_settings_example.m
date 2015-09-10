% Advanced settings file for rtdgui
% EXAMPLE (copy and modify according to requirements)
%
% This file is a script that generates a 'customOpt' struct that is
% equivalent to the 'options' struct that can be passed as an optional
% parameter to rtdTool or rtdPlots.
%
% Allowed entries with default values/behaviors in parentheses 
% (1/0 for true/false settings):
%
% fileList - specify trace file list in getFiles format (prompt user)
% constants - struct to overwrite any cascadeConstants settings (n/a)
%
% scaleAcc - scale acceptor by fixed factor (0)
% skipCriteria - skip filtering by autotrace (0)
% simplePostSync - use non-iterative, thresholded post-synchronization (0)
% idlTraces - idealize traces with SKM (1)
% skipStateFilter - don't divide traces into productive/non-productive (0)
% reidlTraces - reidealize productive/non-productive trace files individually (0)
% stateOcc - generate and display state occupancy over time plots (1)
%
% minFrames - minimum length of events for post-synchronization (35)
% preFrames - frames to keep before point of post-synchronization (5)
% totalFrames - total number of frames in post-synchronized file (500)
% simpleThresh - FRET threshold for non-iterative post-synchronization (0.2)
%
% kinModel - path to QuB model for SKM ('\tRNA selection\2014_04_18 EColi.qmf')
% prodState - state number in kinModel used to identify productive events (4)
% prodDwell - minimum dwell time in milliseconds for productive events (100)
%
% autoPostfix - file name postfix for auto-filtered traces ('_auto')
% selectPostfix - file name postfix for productive traces ('_sel')
% rejectPostfix - file name postfix for non-productive traces ('_rej')
%
% skmOpt - struct to overwrite any SKM default settings (n/a)

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% EXAMPLE SETTINGS

% Don't separate traces into productive and non-productive subpopulations
% (this overwrites the GUI setting)
customOpt.skipStateFilter = 1;              % default = 0

% Display more frames before the point of post-synchronization
customOpt.preFrames = 10;                   % default = 5

% Use a different postfix for traces selected by autotrace criteria
customOpt.autoPostfix = '_filtered';        % default = '_auto'

% Increase the number of iterations for SKM idealization
customOpt.skmOpt.maxItr = 500;              % default = 100

