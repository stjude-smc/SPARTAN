function demux_spots(files, PLT)
% Select molecules with positions that are within printed spots of the
% field of view in multiplexed surface patterning experiments.
% This implementation uses edge detection and circle fitting for spot selection.
% 
% Parameters:
%   files: A string or cell array of input file names.
%   PLT: Plotting behavior
%       - false: No plotting
%       - [] (default): Plot and save using input file basename + '_demux.png'
%       - string: Save plot to the specified file name
%       - cell array: Save plots to corresponding file names for each input file.

% Parameters
suffix = {'A', 'B', 'C', 'D'}; % Add to split file names

% Default value for PLT
if nargin < 2
    PLT = [];
end

% Check input arguments
if nargin < 1
    filter = {'*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
        '*.traces','Binary Traces Files (*.traces)' ; ...
         '*.*','All Files (*.*)'};
    files = getFiles(filter);
elseif ischar(files)
    files = {files};
elseif ~iscell(files)
    error('Input must be a string or cell array of strings');
end

% Validate PLT if it is a cell array
if iscell(PLT) && numel(PLT) ~= numel(files)
    error('PLT must have the same length as the number files if provided as a cell array');
end
if ischar(PLT) && numel(files) > 1
    error('PLT must have the same length as the number of files');
end


for i = 1:numel(files)

    file = files{i};
    output = cell(size(suffix));
    [p, f, e] = fileparts(file);
    % Generate filenames for output files, if not provided
    for j = 1:numel(suffix)
        output{j} = fullfile(p, [f, '_', suffix{j}, e]);
    end
    
    if iscell(PLT)
        plt = PLT{i};
    else
        plt = PLT;
    end
    demux_spots_internal(file, output, plt);
end
end
