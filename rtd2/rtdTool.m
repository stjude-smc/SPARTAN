function [plotWindow] = rtdTool(varargin)
% RTDTOOL Analysis tool for real-time delivery (RTD) experiments
%
%
% RTDTOOL()
% Prompts user to select trace files and processes them with default settings.
%
%
% RTDTOOL(options)
% Use the struct 'options' for user-defined settings. Allowed entries with 
% default values/behaviors in parentheses (1/0 for true/false settings):
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
%
%
% plotWindow = RTDTOOL()
% plotWindow = RTDTOOL(options)
% Returns a handle, plotWindow, to the generated plot window.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


plotWindow = 0;

fprintf('rdTool started.\n\n');
% Initialize timer;
tic;

% Get default settings.
opt = defaultOptions();

% Overwrite defaults with any user-specified settings.
if nargin
    % Recursively merge structs with function from file exchange.
    opt = mergestruct(opt,varargin{1});
end

% Verify the model file exists.
if ~exist(opt.kinModel,'file') && ~opt.plotOnly,
    % Try the current directory first.
    [~,f,e] = fileparts(opt.kinModel);
    opt.kinModel = [f e];
    
    % If not in the current directory, ask the user to find it.
    if ~exist(opt.kinModel,'file'),
        [f,p] = uigetfile( opt.kinModel, 'Load tRNA selection kinetic model' );
        if f==0, return; end  %user hit cancel
        opt.kinModel = fullfile(p,f);
    end
end

% Can't calculate state occupancy without idealization.
if ~opt.idlTraces
    opt.stateOcc = 0;
end

% Prompt user for trace files if not specified.
if ~isfield(opt, 'fileList')
    filter = {'*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
              '*.traces','Binary Traces Files (*.traces)'; ...
              '*.*','All Files (*.*)'};
    prompt = 'Select trace files. Hit cancel when finished.';
    opt.fileList = getFiles(filter, prompt);
end

% Exit if no files selected.
if isempty(opt.fileList)
    fprintf('No trace file selected. Exiting.\n\n');
    return;
end

% Run plot routine if plot-only mode is active, then exit.
if opt.plotOnly
    fprintf('Plot-only mode.\n');
    opt.allFiles = opt.fileList;
    
    % Check if all specified trace files have a corresponding stateOcc file
    counter = 0;
    for i=1:length(opt.allFiles)    
        [path,name,~] = fileparts(opt.allFiles{i});
        opt.stateOccFiles{i} = fullfile(path, [name '.qub_stateOcc.txt']);
        if exist(opt.stateOccFiles{i},'file')
            counter = counter + 1;
        end
    end
    if counter == length(opt.allFiles)
        opt.stateOcc = 1;
    else
        opt.stateOcc = 0;
    end
    
    % Display plots and exit.
    plotWindow = displayPlots(opt);
    return;
end

% Scale acceptor if specified.
if opt.scaleAcc
    for i=1:length(opt.fileList)
        fprintf('%s',['Scaling acceptors: ' opt.fileList{i}]);
        fprintf('\n\n');
        
        [path,name,~] = fileparts(opt.fileList{i});
        outFile = fullfile(path, [name '_corr(' num2str(opt.scaleAcc) ').traces']);
        scaleacceptor(opt.fileList{i}, opt.scaleAcc, outFile);
        
        % Assemble new file list for further processing.
        opt.fileList{i} = outFile;
    end
end

% Select and save traces.
if ~opt.skipCriteria
    for i=1:length(opt.fileList)
        fprintf('%s',['Filtering traces: ' opt.fileList{i}]);
        fprintf('\n\n');
        
        [path,name,~] = fileparts(opt.fileList{i});
        filterOpt.outFilename = fullfile(path, [name opt.autoPostfix '.traces']);
        filterOpt.showWaitbar = 0;
        
        % loadPickSaveTraces() is the core function used by autotrace
        loadPickSaveTraces({opt.fileList{i}}, opt.constants.defaultAutotraceCriteria, filterOpt);
        opt.fileList{i} = filterOpt.outFilename;
    end
end

% Post-synchronize traces.
if ~opt.simplePostSync % Iterative post-synchronization with 2-state model.
    fprintf('Post-synchronizing traces using iterative 2-state idealization method.\n');
    postSyncTraces(opt.minFrames, opt.preFrames, opt.totalFrames, opt.fileList);
    
    % Assemble new file list for further processing.
    for i=1:length(opt.fileList)
        [path,name,~] = fileparts(opt.fileList{i});
        opt.fileList{i} = fullfile(path, [name '_postSync.traces']);
    end
else % Simple threshold-based post-synchronization.
    fprintf('Post-synchronizing traces using non-iterative thresholding method.\n');
    simplePostSync(opt.preFrames, opt.totalFrames, opt.simpleThresh, opt.fileList)
    
    % Assemble new file list for further processing.
    for i=1:length(opt.fileList)
        [path,name,~] = fileparts(opt.fileList{i});
        opt.fileList{i} = fullfile(path, [name '_simplePostSync.traces']);
    end
end

% Perform idealization and related operations.
if opt.idlTraces
    % Load QuB model and prepare for SKM.
    kinModel = QubModel(opt.kinModel);
    kinModel.fixMu    = true( kinModel.nStates,1 );
    kinModel.fixSigma = true( kinModel.nStates,1 );
    fretModel = [to_col(kinModel.mu) to_col(kinModel.sigma)];
    
    for i=1:length(opt.fileList)
        % Idealize traces using SKM.
        fprintf('\n');
        fprintf('%s',['Idealizing traces: ' opt.fileList{i}]);
        fprintf('\n');
        
        currentTraces = loadTraces(opt.fileList{i});
        idl = skm(currentTraces.fret, currentTraces.sampling, kinModel, opt.skmOpt);
        [dwt,offsets] = idlToDwt(idl);
        
        % Save idealizations.
        [path,name,~] = fileparts(opt.fileList{i});
        dwtFile = fullfile(path, [name '.qub.dwt']);
        saveDWT(dwtFile, dwt, offsets, fretModel, currentTraces.sampling);
        
        % Separate traces into "successful" and "unsuccessful" events.
        if ~opt.skipStateFilter
            fprintf('\nDividing traces into productive and non-productive events.\n');
            
            prodFrames = opt.prodDwell/currentTraces.sampling;
            selFile = fullfile(path, [name opt.selectPostfix '.traces']);
            rejFile = fullfile(path, [name opt.rejectPostfix '.traces']);
            
            % Separate traces.
            stateMin(opt.fileList{i}, dwtFile, opt.prodState, prodFrames,selFile,rejFile);
            
            if opt.reidlTraces
                % Re-idealize selected traces using SKM.
                fprintf('\n');
                disp(['Re-idealizing traces: ' selFile]);
                
                currentTraces = loadTraces(selFile);
                idl = skm(currentTraces.fret, currentTraces.sampling, kinModel, opt.skmOpt);
                [dwt,offsets] = idlToDwt(idl);
                
                % Save idealizations.
                dwtFile = fullfile(path, [name opt.selectPostfix '.qub.dwt']);
                saveDWT(dwtFile, dwt, offsets, fretModel, currentTraces.sampling);

                % Re-idealize rejected traces using SKM.
                fprintf('\n');
                disp(['Re-idealizing traces: ' rejFile]);
                
                currentTraces = loadTraces(rejFile);
                idl = skm(currentTraces.fret, currentTraces.sampling, kinModel, opt.skmOpt);
                [dwt,offsets] = idlToDwt(idl);
                
                % Save idealizations.
                dwtFile = fullfile(path, [name opt.rejectPostfix '.qub.dwt']);
                saveDWT(dwtFile, dwt, offsets, fretModel, currentTraces.sampling);
            end
            
            % Assemble new file list for further processing.
            opt.allFiles{(i-1)*3+1} = opt.fileList{i};
            opt.allFiles{(i-1)*3+2} = selFile;
            opt.allFiles{(i-1)*3+3} = rejFile;
        else
            % Copy original file list for further processing.
            opt.allFiles = opt.fileList;  
        end               
    end    
else
    % Copy original file list for further processing.
    opt.allFiles = opt.fileList;  
end

% Generate state occupancy plots.
if opt.stateOcc
    fprintf('\nGenerating state occupancy plots.\n');
    for i=1:length(opt.allFiles)
        [path,name,~] = fileparts(opt.allFiles{i});
        dwtFile = fullfile(path, [name '.qub.dwt']);
        stateOccupancy(dwtFile,opt.totalFrames);
        opt.stateOccFiles{i} = fullfile(path, [name '.qub_stateOcc.txt']);
        fprintf('%s',['Saving to ' opt.stateOccFiles{i}]);
        fprintf('\n');
    end
    fprintf('\n');
end

% Display plots.
plotWindow = displayPlots(opt);

% Display completion time.
fprintf('%s',['rtdTool completed in ' num2str(toc) ' seconds.']);
fprintf('\n');

end % function rtdTool()


function [plotWindow] = displayPlots(opt)
    
% Create plot window with subplots.
nCol = length(opt.allFiles);
if opt.stateOcc
    nRow = 4;
elseif ~opt.idlTraces
    nRow = 2;
else
    nRow = 3;
end
plotWindow = figure;
for i=1:nRow
    for j=1:nCol
        figAxesCell{j,i} = subplot(nRow,nCol,(i-1)*nCol+j, 'Parent',plotWindow);
    end
end
opt.constants.defaultMakeplotsOptions.targetAxes = figAxesCell;

% Create plot titles
plotTitles = cell( length(opt.allFiles),1 );

if opt.skipStateFilter
    for i=1:length(opt.allFiles)
        plotTitles{i} = num2str(i);
    end
else
    for i=1:length(opt.allFiles)/3
        plotTitles{(i-1)*3 + 1} = ['File ' num2str(i) ' all'];
        plotTitles{(i-1)*3 + 2} = ['File ' num2str(i) ' prod.'];
        plotTitles{(i-1)*3 + 3} = ['File ' num2str(i) ' non-prod.'];
    end
end


% Display contour plots, state histograms, and transition density plots.
makeplots(opt.allFiles,plotTitles,opt.constants.defaultMakeplotsOptions);

% Display state occupancy plots if specified.
if opt.stateOcc
    for i=1:nCol
        ax = subplot(nRow,nCol,(nRow-1)*nCol+i, 'Parent',plotWindow);
        stateOcc = dlmread(opt.stateOccFiles{i});
        
        hold(ax,'on');
        for j=2:size(stateOcc,1)
            if (j-1) <= length(opt.stateColors)
                color = opt.stateColors{j-1};
            else
                color = '';
            end
            plot(ax, 100*stateOcc(j,:),color,'LineWidth',1.5);
        end
        if i==1 % make axis labels and formatting consistent with makeplots
            xlim(ax, [0 opt.constants.defaultMakeplotsOptions.contour_length]);
            ylim(ax, [-5 105]);
            xlabel(ax, 'Time (frames)');
            ylabel(ax, 'Occupancy (%)');
        end
        set(ax,'YGrid','on');
        box(ax,'on');
        hold(ax,'off');
    end
end

% Make plots scale simultaneously.
for i=1:nRow
    for j=1:nCol
        figAxes(j) = figAxesCell{j,i};
    end
    linkaxes(figAxes,'xy');
end

% Display plot legend (mapped by number).
fprintf('\nPlot legend:\n');
for i=1:length(opt.allFiles)
    fprintf('%s',[plotTitles{i} ' - ' opt.allFiles{i}]);
    fprintf('\n');
end
fprintf('\n');

end % function displayPlots()


function [defaultOpt] = defaultOptions()

% Specify default settings.
defaultOpt.constants = cascadeConstants();  % global default settings

% rtd-specific filter criteria (type 'traceStat' for definitions).
defaultOpt.constants.defaultAutotraceCriteria.max_firstFret = 0.2;

% Plot settings.
defaultOpt.constants.defaultMakeplotsOptions.truncate_tdplot = true;
defaultOpt.constants.defaultMakeplotsOptions.tdp_max = 0.075;
defaultOpt.constants.defaultMakeplotsOptions.contour_length = 50;
defaultOpt.constants.defaultMakeplotsOptions.normalize = 'total time';
defaultOpt.stateColors = {'k','g','b','r','c','m','y'}; % colors for state occupancy vs. time plot
defaultOpt.constants.defaultMakeplotsOptions.colors = [[0 1 0]; [0 0 1]; ...
    [1 0 0]; [0 1 1]; [1 0 1]; [1 1 0]]; % use the same color scheme for state occ. vs. FRET plot

% General settings.
defaultOpt.scaleAcc = 0;                    % don't scale acceptor
defaultOpt.skipCriteria = 0;                % don't skip autotrace
defaultOpt.simplePostSync = 0;              % use iterative post-sync
defaultOpt.idlTraces = 1;                   % idealize traces
defaultOpt.skipStateFilter = 0;             % filter for productive traces
defaultOpt.reidlTraces = 0;                 % reidealize each file
defaultOpt.stateOcc = 1;                    % state occupancy plots

defaultOpt.minFrames = 35;                  % min. event length for post-sync
defaultOpt.preFrames = 5;                   % frames before point of post-sync
defaultOpt.totalFrames = 500;               % total post-sync'ed trace length
defaultOpt.simpleThresh = 0.2;              % threshold for simple post-sync

defaultOpt.kinModel = '2014_04_18 EColi.qmf';        % QuB model for SKM
defaultOpt.prodState = 4;                   % productive state number
defaultOpt.prodDwell = 120;                 % productive dwell (ms)

defaultOpt.autoPostfix = '_auto';           % autotrace file name postfix
defaultOpt.selectPostfix = '_sel';          % productive file name postfix
defaultOpt.rejectPostfix = '_rej';          % non-prod. file name postfix

% SKM parameters
defaultOpt.skmOpt.maxItr = 100;
defaultOpt.skmOpt.convLL = 0.01;
defaultOpt.skmOpt.zeroEnd = 1;
defaultOpt.skmOpt.seperately = 1;
defaultOpt.skmOpt.quiet = 1;
defaultOpt.skmOpt.fixKinetics = 1;

% Used by rtdPlots() - user normally shouldn't need to specify this
defaultOpt.plotOnly = 0;                    % only plot, don't process


end % function defaultOptions()
