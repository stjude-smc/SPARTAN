function [outFilename,picks,allStats,dataAll]=loadPickSaveTraces( varargin )
% loadPickSaveTraces   Select traces passing defined criteria and save to disk.
% 
%   [] = loadPickSaveTraces( FILES, CRITERIA, ... )
%   Loads FILES, selects traces according to selection CRITERIA (pickTraces.m).
%   FILES may be a filename (string), directory (string), and cell array of
%   filenames. If FILES is a cell array, all of the specified filenames are
%   loaded as one large dataset. If a dictory is specified, all FRET data
%   files (.traces) in this directory will be loaded together. Traces that
%   pass the selection CRITERIA are then saved to disk with the file
%   extension "_auto.traces". A log file is also created.
%
%   TODO: optional argument to pass a list of indexes of traces to load OR the
%   output of traceStat so not all trace data has to be reloaded! This will
%   allow processing of very large datasets in autotrace.
%
%   loadPickSaveTraces( ...,'ParameterName',ProptertyValue,... )
%   The following option parameters (specified with the above format):
%    * showWaitbar  - set to 1 to show (0 to not show) progress indicators.
%    * outputFilename
%    * indexes   - instead of loading the traces and calculating stats, use
%                  these picks and only load what is needed. This is a cell
%                  array with one per file.
%    * stats

% TODO:
%  * Allow user to give stats of the data files, so that they aren't
%    calculated a second time in autotrace.


% USER TUNABLE PARAMETERS
options.batchMode = 0;
options.showWaitbar = 1;


%% ---- Process command arguments

if nargin<2,
    error('Must specify files and criteria values');
end
criteria = varargin{2};


% If a directory name is give, find all traces files in the directory.
if ischar(varargin{1}) && isdir(varargin{1}),
    dirName = varargin{1};
    files = dir( [dirName filesep '*.traces'] );

    if numel(files)<1
        disp('No files in this directory!');
        return;
    end

    files = strcat( [dirName filesep], {files.name} );

else
    files = varargin{1};
    
    % If only 1 file given, make into a cell array for consistency.
    if ischar(files), files={files}; end  
    
    % Find a directory common to all files (for logging)
    dirName = commonDir(files);
end

% Verify that the input is a file list
assert( iscell(files), 'Invalid file list' );
nFiles = numel(files);


% Process optional argument list
if nargin>2,
    vopt = varargin{3};
%     vopt = varargin(3:end);
%     assert( mod(numel(vopt),2)==0 & all( cellfun(@ischar,vopt(1:2:end)) ), ...
%             'Invalid format for optional argument list' );
    
    % Merge default options and those specified by the user.
    % The user's options will override any existing defaults.
    options = catstruct( options, vopt );
end


% Verify arguments
if options.batchMode,
    error('Batch mode not yet supported!');
end



%% ---- Load each file, select traces, and add it to output.


% Open file handles for output
if isfield(options,'outFilename'),
    outFilename = options.outFilename;
else
    outFilename = strrep(files{1}, '.traces', '_auto.traces');
    outFilename = strrep(outFilename, '_01_auto.traces', '_auto.traces');
end


% Process each file individually, adding the selected traces to output.
if options.showWaitbar,
    wbh = waitbar(0,'Selecting and saving traces according to criteria...');
end


% Ideally, space should be pre-allocated. TODO.
allStats = struct( [] );
picks = [];
nTracesPerFile = zeros(nFiles,1);
channelNames = {};
        
for i=1:nFiles,
    % Load traces file data. Load only the useful subset if indexes are
    % given to save memory.
    if isfield(options,'indexes'),
        data = loadTraces( files{i}, options.indexes{i} );
    else
        data = loadTraces( files{i} );
    end
    
    % Verify that all files have the same imaging geometry/channel definitions.
    if i==1,
        channelNames = data.channelNames;
    else
        assert( all(strcmp(sort(channelNames),sort(data.channelNames))), 'Traces data format mismatch (different geometry?)' );
    end
  
    % Unless provided, calc trace properties and select those that meet criteria.
    if isfield( options, 'indexes' ),
        % Only the picked traces are loaded above. Nothing to do.
        %indexes = 1:numel(options.indexes{i});
        %stats = struct([]);
    else
        % Calculate trace statistics
        % FIXME: this should check if stats are being passed!!
        if isfield( options, 'stats' ),
            stats = options.stats( sum(nTracesPerFile)+(1:size(data.donor,1)) );
        else
            stats = traceStat( data );
        end
        
        % Pick traces passing criteria
        indexes = pickTraces( stats, criteria );
        indexes = reshape(indexes,1,numel(indexes));  %insure row vector shape.
        picks = [picks indexes+sum(nTracesPerFile)];  %indexes of picked molecules as if all files were concatinated.
        
        data = data.subset(indexes);
        
        % Save stats info for logging.
        if i==1
            allStats = stats;
        else
            allStats = cat(2, allStats, stats);
        end
    end
    
    % Add selected traces from the current file to the output
    if i==1,
        dataAll = data;
    else
        dataAll.combine(data);
    end
    
    % Update status bar
    if options.showWaitbar,
        waitbar(i/nFiles,wbh);
    end
    
    nTracesPerFile(i) = size(data.donor,1);
end

if isfield(options,'stats'),
    allStats = options.stats;
end


% Save trace data to file
saveTraces(outFilename,dataAll);

if options.showWaitbar,
    waitbar(1,wbh,'Finished!');
end



%% ---- Save log file
[p,n,ext] = fileparts(outFilename);
if isempty(p)
    p = pwd;
end
    
% logFilename = strrep(outFilename,'.traces','.log');
logFilename = fullfile(p, [n '.log']);
fid = fopen(logFilename,'w');

% Save header, with filenames and the number of traces picked.
fprintf(fid,'%s\n\n%s\n',date,'DIRECTORY');
fprintf(fid,'  %s\n\n%s\n',dirName,'FILES');

for i=1:numel(files)
    shortName = files{i}(numel(dirName)+1:end);
    fprintf(fid,' %5d: %s\n',nTracesPerFile(i),shortName);
end

% nPicked = numel(picks);
nPicked = size( data.donor,1 );
nTotalTraces = sum(nTracesPerFile);

fprintf(fid,'\nMolecules Picked:\t%d of %d (%.1f%%)\n\n\n', ...
            nPicked, nTotalTraces, 100*nPicked/nTotalTraces );  

        
% Descriptive statistics about dataset
stats = allStats;
isMolecule      = sum( [stats.snr]>0 );
singleMolecule  = sum( [stats.snr]>0 & [stats.overlap]==0 );
hasFRET         = sum( [stats.snr]>0 & [stats.overlap]==0 & [stats.acclife]>=5 );
other           = nPicked;

fprintf(fid,'PICKING RESULTS\n');
fprintf(fid, '  %20s:  %-5d (%.1f%% of total)\n', 'Donor photobleaches', isMolecule,     100*isMolecule/nTotalTraces);
fprintf(fid, '  %20s:  %-5d (%.1f%% of above)\n', 'Single donor',        singleMolecule, 100*singleMolecule/isMolecule);
fprintf(fid, '  %20s:  %-5d (%.1f%% of above)\n', 'Have FRET',           hasFRET,        100*hasFRET/singleMolecule);
fprintf(fid, '  %20s:  %-5d (%.1f%% of above)\n', 'Pass all criteria',   other,          100*other/hasFRET);
fprintf(fid, '\n\n');


% Save picking criteria used
fprintf(fid,'PICKING CRITERIA\n');

names = fieldnames(  criteria );
vals  = struct2cell( criteria );

for i=1:numel(names),
    if isempty( vals{i} ), continue; end  %skip unchecked criteria
    fprintf(fid, '  %22s:  %.2f\n', names{i}, vals{i});
end

% Save values of all other constants used
fprintf(fid, '\n\nCONSTANTS\n');

constants = cascadeConstants;
names = fieldnames(  constants );
vals  = struct2cell( constants );

for i=1:numel(names),
    if isstruct( vals{i} ) || numel( vals{i} )>1, continue; end
    
    fprintf(fid, '  %22s:  %.2f\n', names{i}, vals{i});
end


fprintf(fid,'\n\n');
fclose(fid);


if options.showWaitbar,
    close(wbh);
    drawnow;
end





end




%%
function dirName = commonDir( files )
% Find the common directory containing all give FILES.

% Take out only path names from files
nFiles = numel(files);
for i=1:nFiles,
    files{i} = fileparts(files{i});
end

% Concatinate all pathnames into a single string matrix.
files = strvcat(files);

% Find any differences
diffs = zeros(1,size(files,2));
for i=1:nFiles
    diffs = diffs | files(i,:)~=files(1,:);
end

lastDiff = find(diffs);
if isempty(lastDiff), lastDiff = size(files,2); end

% Find directory name by going back to the last path seperating character.
dirName = fileparts( files(1,1:lastDiff) );


end

