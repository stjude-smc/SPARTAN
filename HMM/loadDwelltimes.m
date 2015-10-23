function [dwellc,sampling,model] = loadDwelltimes( filename, options )
% loadDwelltimes  Combined list of dwell-times for each state.
%
%   [DWELLS,SAMPLING,MODEL] = loadDwelltimes(FILENAME) creates the cell array
%   DWELLS with a list of all dwell times (in seconds) from each state. 
%   SAMPLING is the time resolution in ms and MODEL is the FRET model.
%
%   [...] = loadDwelltimes(...,'removeBlinks') will remove dwells in the 
%   zero-FRET state (state 1) associated with blinking.
%
%   See also: removeBlinks, dwellhist, lifetime_exp.

%   Copyright 2007-2015 Cornell University All Rights Reserved.



%% Prompt user for file names if not given.
if nargin<1,
    filename = getFile('*.dwt','Choose a dwell-time file');
    if isempty(filename), return; end  %user hit cancel.
end


% If .traces files are given, quietly look for the corresponding .qub.dwt
[p,f,e] = fileparts(filename);
    
if ~strcmp(e,'.dwt'),
    filename = fullfile(p,[f '.qub.dwt']);
    if ~exist(filename, 'file'),
        filename = fullfile(p,[f '.dwt']);
    end
    if ~exist(filename,'file'),
        error('Input file is not a .dwt and no associated .dwt file could be found');
    end
end



%% Load the dwell times and remove dark-state dwells.
[dwells,sampling,offsets,model] = loadDWT(filename);
assert( numel(dwells)>0, 'Empty or invalid dwell-time file' );

nStates = numel(model)/2;

if nargin>=2 && strcmpi(options,'removeBlinks'),
    dwells = removeBlinks(dwells,offsets);
end

% Merge dwells from all traces
dwells = vertcat(dwells{:});
states = dwells(:,1);
times = dwells(:,2);

% Build list of dwell times in each class
dwellc = cell(nStates, 1);

for i=1:nStates,
    dwellc{i} = times(states==i)*sampling/1000;
end



