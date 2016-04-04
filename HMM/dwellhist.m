function varargout = dwellhist( dwtfilename, inputParams )
%dwellhist  Dwell-time histograms
% 
%   dwellhist(FILES) displays log-scale dwell-time histograms for each
%   .dwt file in the cell array FILES (must all have same time resolution).
%
%   HIST = dwellhist(FILES) returns histograms in the columns of HIST as:
%     [X_axis, State1/File1, State2/File1, ..., State1/File2, State2/File2, ...]
%   The X-axis is in seconds.
%
%   dwellhist() prompts the user for a list of files.
%
%   dwellhist(...,PARAMS) give optional parameters in struct PARAMS:
%
%      'removeBlinks': Remove dwells in dark states (class 1). Default=true.
%                      Dwells broken up by such blinks are merged, with the
%                      time during the blink added to surrounding dwells.
%
%      'logX':         Use a log-scale times axis to aid visualization of 
%                      multi- exponential distributions (default=true).
%                      See Sigworth and Sine (1987), Biophys J 50, p. 1047-1054.
%
%      'dx':           Log-scale time axis bin size. default=0.2.
%
%      'normalize':    log-scale histogram normalization method:
%                      'off'   - raw dwell counts, no normalization.
%                      'state' - each state; sum of each histogram=1. (default)
%                      'file'  - all states; sum of all histograms per file=1.
%                      'time'  - dwell counts per second of observation time
%
%   See also: lifetime_exp, loadDwelltimes, removeBlinks.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


%% ---- USER TUNABLE PARAMETERS ----

params.logX = true;
params.dx = 0.2;  %log-scale bin width (0.1=25%, 0.2=60%, 0.5=3-fold, 1=10-fold)
params.removeBlinks = true;
params.normalize = 'state';

% Merge options, giving the user's options precedence.
if nargin>1,
    params = mergestruct( params, inputParams );
end

% Check parameters
assert( all(ismember(params.normalize,{'off','state','file','time'})), ...
        'Invalid normalization option' );



%% Prompt user for file names if not given.

% Get list of .dwt files to load
if nargin<1,
    dwtfilename = getFiles('*.dwt','Choose dwell-time files');
end
dwtfilename = findDwt(dwtfilename);
names = trimtitles(dwtfilename);
if numel(dwtfilename)==0, return; end


% If the histogram output is requested, don't plot anything.
if nargout>0,
    [dwellaxis,histograms,meanTime] = dwellhist2(dwtfilename,params);
    dwellhist = [to_col(dwellaxis) horzcat(histograms{:})];
    
    output = {dwellhist,names,meanTime};
    [varargout{1:nargout}] = output{1:nargout};
    
% Otherwise, plot dwell times
else
    % Save the histogram data in the figure for access by callback functions.
    hFig = figure;
    handles.hFig = hFig;
    handles.names = names;
    handles.params = params;
    handles.dwtfilename = dwtfilename;
    guidata(hFig,handles);

    % Add a control at the bottom of the GUI for saving the histograms to file.
    uicontrol( 'Style','pushbutton', 'String','Save...', ...
               'Position',[15 15 75 30], 'Callback',@saveDwelltimes, ...
               'Parent',hFig );
           
%     uicontrol( 'Style','pushbutton', 'String','Settings...', ...
%                'Position',[115 15 75 30], 'Callback',@dwellhist_settings, ...
%                'Parent',hFig );
           
    uicontrol( 'Style','pushbutton', 'String','Replot...', ...
               'Position',[215 15 75 30], 'Callback',@dwellhist_plot, ...
               'Parent',hFig );
       
    dwellhist_plot(hFig);
end



end %function dwellhist





%% ------ Load dwell-times from .dwt files
function [dwellaxis,histograms,meanTime] = dwellhist2(dwtfilename,params)

nFiles = numel(dwtfilename);
dwells  = cell(nFiles,1);  %consolidated list of dwell times in each state
sampling = zeros(nFiles,1);

for i=1:nFiles,
    % Load dwell-times a list per state, concatinating dwells from all traces,
    % ignoring the zero-FRET state if applicable.
    if params.removeBlinks,
        [dwells{i},sampling(i)] = loadDwelltimes( dwtfilename{i}, 'removeBlinks' );
        dwells{i} = dwells{i}(2:end);
    else
        [dwells{i},sampling(i)] = loadDwelltimes( dwtfilename{i} );
    end
end

% Verify all input idealizations are roughly consistent.
if ~all(sampling==sampling(1)),
    warning('Mismatched time resolution of input files');
end
sampling = sampling(1)/1000;  %convert to seconds.

nStates = cellfun(@numel, dwells);
if ~all(nStates==nStates(1)),
    warning('Idealization models have a different number of states');
end
nStates = max(nStates);


% Get dwell time limits for setting axes limits later.
meanTime = zeros(nFiles,nStates);  %mean dwell-time per state/file.
maxTime = 0;  %longest dwell in seconds
totalTime = zeros(nFiles,nStates);

for i=1:nFiles,
    dwellc = dwells{i};
    maxTime = max( maxTime, max(vertcat(dwellc{:})) );
    meanTime(i,:) = cellfun(@mean, dwellc)';
    totalTime(i,:) = cellfun(@sum, dwellc)';
end



%% Calculate dwell time bins (EDGES)
if ~params.logX,
    % Linear X-axis in seconds.
    dwellaxis = 0:sampling:maxTime;
else
    % Create a log time axis with a fixed number of bins.
    % histcounts uses bin edges of E(k) <= X(i) < E(k+1).
    dwellaxis = log10(sampling):params.dx:log10(maxTime*3);
    
    % Force the bins edges to be exact intervals of the time resolution.
    % The histogram will better sample the discrete nature of the data.
    maxFrames = ceil(maxTime*3/sampling);
    fullaxis = log10( (1:maxFrames)*sampling )';
    dwellaxis = unique( nearestBin(dwellaxis, fullaxis) );
%     dwellaxis = unique( floor(fullaxis/sampling)*sampling );
    
    % Normalization factor to account for varying-sized bins.
    dlx = dwellaxis(2:end) - dwellaxis(1:end-1);
    dlx = [dlx dlx(end)];
end


%% Calculate histograms
histograms = cell(nFiles,nStates);

for file=1:nFiles,
    ndwells = cellfun(@numel,dwells{file});
            
    for state=1:nStates,
        % Small constant ensures dwells fall in the correct histogram bin.
        dwellc = dwells{file}{state} +sampling/10;
        
        % Make linear-scale survival plot.
        if ~params.logX,
            counts = histc( dwellc, dwellaxis );
            histdata = sum(counts) - cumsum(counts);
            histdata = histdata/histdata(1);
        
        % Make log-scale Sine-Sigworth plot (linear ordinate).
        else
            counts = histc( log10(dwellc)', dwellaxis );
            histdata = counts./dlx;  %normalize by log-space bin size
            histdata = histdata/sum(histdata);  %normalize to 1
            
            switch params.normalize
                case 'off'  %raw dwell counts
                    histdata = histdata*ndwells(state);
                case 'state'  %fraction of counts in each bin for this state
                    histdata = 100*histdata;
                case 'file'  %fraction of counts in each bin across entire file
                    histdata = 100*histdata *ndwells(state)/sum(ndwells);
                case 'time'  %fraction of dwells per total observation time
                    histdata = histdata*ndwells(state)/sum(totalTime(file,:));
            end
        end
        
        histograms{file,state} = to_col(histdata);
    end
end


% Combine histograms into a matrix for saving.
if params.logX,
    dwellaxis = 10.^dwellaxis;
end



end %function dwellhist2



function [newVal,idx] = nearestBin( values, bins )
% For each VALUE, find the BIN with the closest value.

newVal = zeros( size(values) );
idx = zeros( size(values) );

for i=1:numel(values),
    [~,idx(i)] = min( abs(bins-values(i)) );
    newVal(i) = bins(idx(i));
end

end





%% Display the histograms
function dwellhist_plot(hObject,~,~)
% Actually plots the dwell-time histograms.


% Get GUI data from initial call
handles = guidata(hObject);
names   = handles.names;
params  = handles.params;

set(handles.hFig,'pointer','watch'); drawnow;


% Recalculate histograms directly from file.
[dwellaxis,histograms] = dwellhist2(handles.dwtfilename,params);

handles.dwellhist = [to_col(dwellaxis) horzcat(histograms{:})];
guidata(hObject,handles);  %update for later calls to saveDwelltimes()


% Find a good zoom axis range for viewing all of the histograms.
% if params.logX,
    xmax = dwellaxis(end);
% else
%     xmax = 4* max(meanTime(:));
% end
h = [histograms{:}];
ymax = max(h(:));


% Choose ordinate label based on normalization
if params.logX
    switch params.normalize
        case 'off'
            ordinate = 'Counts';
        case 'state'
            ordinate = 'Counts (%)';
        case 'file'
            ordinate = 'Counts (% of file)';
        case 'time'
            ordinate = 'Counts s^{-1}';
    end
end


% Display survival plots, one state per panel.
[~,nStates] = size(histograms);
ax = zeros(nStates,1);

for state=1:nStates,
    ax(state) = subplot( nStates, 1, state, 'Parent',handles.hFig );

    if params.logX,
        semilogx( ax(state), dwellaxis, [histograms{:,state}], '-', 'LineWidth',2 );
        ylabel(ax(state), ordinate);
    else
        plot( ax(state), dwellaxis, [histograms{:,state}], '-', 'LineWidth',2 );
        ylabel(ax(state), 'Dwell Survival (%)');
    end
    if state==nStates,
        xlabel(ax(state), 'Time (s)');
    end
    
    xlim( ax(state), [dwellaxis(1) xmax] );
    ylim( ax(state), [0 ymax] );
    title(ax(state), sprintf('State %d',state) );
end

legend(ax(end), names);
linkaxes(ax,'xy');


set(handles.hFig,'pointer','arrow'); drawnow;

end %function dwellhist_plot




%% ------ Save results to file for plotting in Origin
function saveDwelltimes(hObject,~,~)
% Callback function for the "save histograms" button in the histogram figure.

% Get histogram data saved in the figure object.
handles = guidata(hObject);
dwellhist = handles.dwellhist;
names = handles.names;

nFiles = numel(names);
nStates = (size(dwellhist,2)-1)/nFiles;

% Ask the user for an output filename.
[f,p] = uiputfile('*.txt','Select output filename','dwellhist.txt');
if f==0, return; end  %user hit cancel.
outFilename = fullfile(p,f);


% Output header lines
fid = fopen(outFilename,'w');
fprintf(fid,'Time (s)');

for state=1:nStates,
    for i=1:nFiles
        fprintf(fid,'\tState%d %s',state,names{i});
    end
end
fprintf(fid,'\n');
fclose(fid);


% Output histogram data
dlmwrite(outFilename, dwellhist, 'delimiter','\t', '-append');


end


