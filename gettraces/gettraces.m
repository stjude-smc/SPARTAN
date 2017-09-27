function [stkData,peaks] = gettraces(varargin)
% GETTRACES  Extract smFluorescence traces from movies
%
%    GETTRACES() launches the graphical user interface.
%
%    [STK,PEAKS] = GETTRACES( FILENAME, PARAMS )
%    Loads a movie from FILENAME, finds fluorescence peaks, and returns
%    their locations (PEAKS) as a Nx2 matrix (x,y).  IMG is the image used
%    for selecting intensity peaks (sum of aligned fluorescence channels).
%    STK is the stkData structure containing the internal state.
%
%    [STK,PEAKS] = gettraces( FILENAME, PARAMS, OUTFILE )
%    As above, but also extracts traces from the PEAKS locations and saves
%    the resulting traces to OUTFILE.
%
%    STK can also be passed instead of FILENAME if the movie has already
%    been loaded by gettraces. This saves on load time.
%
%    GETTRACES( DIRECTORY, PARAMS )
%    This is "batch mode".  For each filename in DIRECTORY, the movie is loaded,
%    processed, and a .rawtraces file is saved automatically. Additional
%    (optional) fields in PARAMS include option to check all child folders
%    (recursive) and option to skip processed movies (skipExisting).
%
%
%  PARAMS is a struct with any of the following options:
% 
% - don_thresh:     total (D+A) intensity threshold for pick selection
%                   (if 0, threshold is automatically selected).
% 
% - overlap_thresh: peaks are rejected if they are closer than this radius.
% 
% - nPixelsIntegrated: number of pixels proximal to each peak to sum
%      to produce fluorescence traces. Higher values capture more
%      intensity, but also capture more noise...
% 
% - geometry: imaging geometry can be single-channel/full-chip (1), 
%             dual-channel/left-right (2), or quad-channel (3/4).
%             Default: dual-channel (2).
% 
% - crosstalk: donor-to-acceptor channel fluorescence channel
%              bleedthrough (crosstalk) as a fraction. For correction.
% 
% - photonConversion: fluorescence ADU/photon conversion factor.
% 
% - fieldNames: cell array of assignments of fluorescence field names
%                 (e.g., donor, acceptor, factor). If one is left blank
%                 that field is disregarded. OPTIONAL

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% NOTE: because of the structure of this file, all variables defined in this
% initial section are global to the other sub-functions in the file. This way,
% common parameters (params, stkData, constants) do not have to be passed
% around. FIXME: avoid using common names where possible to pervent confusion,
% including image_t, peaks, and picks.


% If calling directly from command line, launch the GUI.
if nargin==0,
    gettraces_gui;
    return;
end


%------ Load parameter values (always second parameter!)
% Set default parameter values and merge with user-specified values.
% The user's options will override any existing defaults.
constants = cascadeConstants;
params = constants.gettracesDefaultParams;

if nargin>=2,
    userParams = varargin{2};
    % FIXME: verify all required parameter settings are given!    
    params = mergestruct( params, userParams );
end


%------ If a structure is specified, load as STK data
if nargin>=1 && isstruct(varargin{1}),
    stkData = varargin{1};

%------ Otherwise, load a list of files (cell array) or single file.
elseif nargin>=1 && iscell(varargin{1}) || (ischar(varargin{1}) && ~isdir(varargin{1})),
    stkData = OpenStk(varargin{1}, params);

%------ If a directory name is given, run in batch mode and then terminate.
elseif nargin>=1 && ischar(varargin{1}) && isdir(varargin{1})
    batchmode(varargin{1}, params);
    return;

else
    error('gettraces: Invalid param 1');
end

% If all we want is the movie data, don't pick yet.
if nargout==1,
    return;
end


%------ Find peaks of intensity corrosponding to isolated single molecules

% Find peak locations from total intensity.
% stkData outputs: total_t, alignStatus, total_peaks, fractionOverlapped, 
%    rejectedPicks, rejectedTotalPicks
stkData = getPeaks(stkData, params);
peaks = stkData.peaks;

% Generate integration windows for later extracting traces.
% Outputs: regions, integrationEfficiency, fractionWinOverlap.
stkData = getIntegrationWindows(stkData, params);

%------ Extract traces using peak locations
if nargin>=3 && ischar(varargin{3}),
    integrateAndSave(stkData, varargin{3}, params)
end
    
return;

end %FUNCTION gettraces




%% --------------------- OPEN STK MOVIE (CALLBACK) --------------------- %
function batchmode(direct,params)

% Get list of files in current directory (option: and all subdirectories)
movieFiles = regexpdir(direct,'^.*\.(tiff?|stk)$',params.recursive);
nFiles = length(movieFiles);

% Wait for 100ms to give sufficient time for polling file sizes in the
% main loop below.
pause(0.1);


% Process each file in the user selected directory.
nTraces  = zeros(nFiles,1); % number of peaks found in file X
existing = zeros(nFiles,1); % true if file was skipped (traces file exists)

for i=1:nFiles,
    
    % Skip if previously processed (.rawtraces file exists)
    stk_fname = movieFiles(i).name;
    [p,name]=fileparts(stk_fname);
    traceFname = fullfile(p, [name '.rawtraces']);
    
    if params.skipExisting && exist(traceFname,'file'),
        if ~params.quiet,
            disp( ['Skipping (already processed): ' stk_fname] );
        end
        existing(i) = 1;
        continue;
    end
    
    % Poll the file to make sure it isn't changing.
    % This could happen when a file is being saved during acquisition.
    d = dir(movieFiles(i).name);
    if movieFiles(i).datenum ~= d(1).datenum,
        disp( ['Skipping (save in process?): ' movieFiles(i).name] );
        existing(i) = 1;
        continue;
    end
    
    % Show waitbar only when new data must be loaded.
    if ~exist('h','var') && ~params.quiet,
        h = waitbar( (i-1)/nFiles,'Extracting traces from movies...');
    end
    
    % Load STK file
    try
        gettraces( movieFiles(i).name, params, traceFname );
    catch e
        disp('Skipping file: corrupted, missing, or not completely saved.');
        disp(movieFiles(i).name);
        disp(e);
        existing(i) = 1;
        continue;
    end
    
    if exist('h','var'),
        waitbar(i/nFiles, h); drawnow;
    end
end

if exist('h','var'),  close(h);  end

% If no new data was loaded, nothing more to do.
if all(existing),  return;  end


% ----- Create log file with results
log_fid = fopen( fullfile(direct,'gettraces.log'), 'wt' );

% Log parameter values used in gettraces
fprintf(log_fid,'GETTRACES PARAMETERS:\n');

output = evalc('disp(params)');
fprintf(log_fid, '%s', output);
%FIXME: structure parameters are not displayed here (alignment!)

% Log list of files processed by gettraces
fprintf(log_fid,'\n%s\n\n%s\n%s\n\n%s\n',date,'DIRECTORY',direct,'FILES');

for i=1:nFiles
    if existing(i),
        fprintf(log_fid, 'SKIP %s\n', movieFiles(i).name);
    else
        fprintf(log_fid, '%.0f %s\n', nTraces(i)/2, movieFiles(i).name);
    end
end


% Clean up
fclose(log_fid);

end %FUNCTION batchmode





