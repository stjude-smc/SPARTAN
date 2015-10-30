function [meanTime,dwellhist,names] = dwellhist( dwtfilename, inputParams )
% dwellhist  Survival plots of state dwell-times
% 
%   [MEANS,HIST] = dwellhist(FILES) creates survival plots of state dwell-times 
%   for each state in each .dwt file in the cell array FILES. The first column 
%   of HIST is the X-axis for all plots. For each state, HIST has a series of
%   consecutivecolumns for each file in FILES so that the columns order is:
%       State1/File1, State1/File2, ... State2/File1, State2/File2, ...
%   MEANS contains mean dwell-times for each state (columns) and file (rows).
%   The mean of an exponential distribution is approximately the time constant.
%
%   dwellhist(FILES) with no output arguments displays the survival plots,
%   with a button to save the histgrams as a text file with column headers.
%
%   [...] = dwellhist(...,PARAMS) specifies optional parameters in the struct
%   PARAMS to control how the histograms are made (true/false values):
%
%      'removeBlinks': Ignore dwells in the dark state, assumed to be state 1.
%                      Dwells broken up by such blinks are merged, with the
%                      time during the blink added to surrounding dwells.
%
%      'logX':         Use a log-scale times axis to aid visualization of 
%                      multi- exponential distributions. See Sigworth and Sine
%                      (1987), Biophys J 50, p. 1047-1054.
%
%   See also: lifetime_exp, loadDwelltimes, removeBlinks.

%   Copyright 2007-2015 Cornell University All Rights Reserved.



%% ---- USER TUNABLE PARAMETERS ----

% Sine-sigworth transformation.
params.logX = true;
nbins = 25;  %only for log-X.

% Remove blinking events (dwells in state 1)
params.removeBlinks = true;  % merge blinks into previous dwell

% Merge options, giving the user's options precedence.
if nargin>1,
    params = mergestruct( params, inputParams );
end



%% Prompt user for file names if not given.
if nargin<1,
    disp('Select DWT files, hit cancel when finished');
    dwtfilename = getFiles('*.dwt','Choose dwell-time files');
end

% if just one filename given, convert to cell array
if ~iscell( dwtfilename ),
    dwtfilename = {dwtfilename};
end

nFiles = numel(dwtfilename);
if nFiles==0, return; end



%% ------ Load dwell-times from .dwt files
names = cell(nFiles,1);
dwts  = cell(nFiles,1);
sampling = zeros(nFiles,1);

for i=1:nFiles,
    % Generate short names for each file for display.
    [~,f] = fileparts(dwtfilename{i});
    names{i} = strrep( strrep(f,'.qub',''), '_',' ' );
    
    % Load dwell-times a list per state, concatinating dwells from all traces,
    % ignoring the zero-FRET state if applicable.
    if params.removeBlinks,
        [dwellc,sampling(i)] = loadDwelltimes( dwtfilename{i}, 'removeBlinks' );
        dwellc = dwellc(2:end);
    else
        [dwellc,sampling(i)] = loadDwelltimes( dwtfilename{i} );
    end
    
    dwts{i} = dwellc;
end

% Verify all input idealizations are roughly consistent.
if ~all(sampling==sampling(1)),
    warning('Mismatched time resolution of input files');
end
sampling = sampling(1)/1000;  %convert to seconds.

nStates = cellfun(@numel, dwts);
if ~all(nStates==nStates(1)),
    warning('Idealization models have a different number of states');
end
nStates = max(nStates);


% Collect dwell statistics for later.
meanTime = zeros(nFiles,nStates);  %mean dwell-time per state/file.
maxTime = 0;  %longest dwell in seconds

for i=1:nFiles,
    dwellc = dwts{i};
    maxTime = max( maxTime, max(vertcat(dwellc{:})) );
    meanTime(i,:) = cellfun(@mean, dwellc)';
end



%% Create dwell-time survival histograms

% Create x-axis bins for histc. Note that these are BIN EDGES!
if ~params.logX,
    % Linear X-axis in seconds.
    dwellaxis = 0:sampling:maxTime;
    dwellhist = zeros( numel(dwellaxis), 1+(nStates*nFiles) );
    dwellhist(:,1) = dwellaxis;
else
    % Create a log millisecond time axis with a fixed number of bins.
    dwellaxis = linspace( log10(sampling), 3*log10(maxTime), nbins );
    
    % Force the bins edges to be exact intervals of the time resolution.
    % The histogram will better sample the discrete nature of the data.
    fullaxis = log10( (1:5000)*sampling )';
    dwellaxis = unique( nearestBin(dwellaxis, fullaxis) );
    
    % Normalization factor to account for varying-sized bins.
    dlx = diff( [dwellaxis(1:end-1) ; dwellaxis(2:end)] );
    dlt = 10.^( [dlx dlx(end)] );
end


for file=1:nFiles,
    for state=1:nStates,
        dwellc = dwts{file}{state};
        
        % Make linear-scale survival plot.
        if ~params.logX,
            counts = histc( dwellc, dwellaxis );
            histdata = sum(counts) - cumsum(counts);
            histdata = histdata/histdata(1);
        
        % Make log-scale Sine-Sigworth plot (linear ordinate, log10 ms X).
        else
            counts = histc( log10(dwellc)', dwellaxis );
            histdata = counts./dlt;  %normalize by log-space bin size
            histdata = histdata/sum(histdata);  %normalize to 1
        end
        
        % Add to the output matrix
        colIdx = ((state-1)*nFiles)+file+1;
        dwellhist( :, colIdx ) = histdata;
    end
end

% Do not display plots or save the output if the histogram matrix is requested.
if nargout>0, return; end

if nFiles>1,
    disp('Mean dwell-times are listed with files across rows and states across columns.');
end



%% Display the histograms

% If there are a ton of datapoints, we can't use the simple colors above.
% Instead, just use a simple blue-to-red gradient
% if nFiles>5,
%     colors = zeros(nFiles,3);
%     interval = (1/(nFiles-1));
%     colors(:,1) = 0:interval:1;
%     colors(:,3) = 1:-interval:0;
% end

% Find a good zoom axis range for viewing all of the histograms.
if ~params.logX,
    xmax = 4* max(meanTime(:));
    xtitle = 'Time (s)';
else
    xmax = dwellaxis(end);
    xtitle = 'Time (Log10 s)';
end

h = dwellhist(:,2:end);
ymax = max( h(:) );

% Save the histogram data in the figure and add a button so the data
% can be saved to file by the user.
hFig = figure;
setappdata(hFig,'dwellhist',dwellhist);
setappdata(hFig,'names',names);

% Display survival plots, one state per panel.
ax = zeros(nStates,1);

for state=1:nStates,
    ax(state) = subplot( nStates, 1, state );
    colIdx = (state-1)*nFiles + (1:nFiles) +1;
    plot( dwellaxis, dwellhist(:,colIdx), '-', 'LineWidth',2 );

    xlim( [dwellaxis(1) xmax] );
    ylim( [0 ymax] );

    if ~params.logX,
        ylabel('Dwell Survival (%)' );
    else
        ylabel('Counts (%)');
    end
    if state==nStates,
        xlabel(xtitle);
    end
    
    title( sprintf('State %d',state) );
end

legend(names);
linkaxes(ax,'xy');

% Add a control at the bottom of the GUI for saving the histograms to file.
uicontrol( 'Style','pushbutton', 'String','Save...', ...
           'Position',[15 15 75 30], 'Callback',@saveDwelltimes, ...
           'Parent',hFig );


end %function dwellhist


function [newVal,idx] = nearestBin( values, bins )
% For each VALUE, find the BIN with the closest value.

newVal = zeros( size(values) );
idx = zeros( size(values) );

for i=1:numel(values),
    [~,idx(i)] = min( abs(bins-values(i)) );
    newVal(i) = bins(idx(i));
end

end



%% ------ Save results to file for plotting in Origin
function saveDwelltimes(hObject,~,~)
% Callback function for the "save histograms" button in the histogram figure.

% Get histogram data saved in the figure object.
hFig = get(hObject,'Parent');
dwellhist = getappdata(hFig,'dwellhist');
names = getappdata(hFig,'names');

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


