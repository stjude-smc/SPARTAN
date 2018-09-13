function varargout = dwellhist(varargin)
%dwellhist  Dwell-time histograms
% 
%   [X,HIST] = dwellhist(FILES) dwell-time histograms for each .dwt file in 
%   the cell array FILES (must all have same time resolution). HIST is a cell
%   array with files in rows and states in columns size=[nFiles,nStates].
%   The X are the bin edges in seconds.
%
%   [X,HIST,FITS] = dwellhist(FILES) also calculates fit lines (FITS) for each 
%   state using the mean dwell time parameters provided in the optional 
%   parameter meanDwellTime (see below). One column in FITS per state.
%   Bins are the same as for histograms.
%
%   dwellhist(...) without output parameters creates a new figure to display
%   the calculated dwell-time histograms.
%
%   dwellhist(FILES,PARAMS) give optional parameters in struct PARAMS:
%
%      'removeBlinks': Remove dwells in dark states (class 1). Default=true.
%                      Dwells broken up by such blinks are merged, with the
%                      time during the blink added to surrounding dwells.
%
%      'logX':         Use a log-scale times axis to aid visualization of 
%                      multi- exponential distributions (default=true).
%                      See Sigworth and Sine (1987), Biophys J 50, p. 1047-1054.
%
%      'dx':           Log-scale time axis bin size. default=0.4.
%
%      'normalize':    log-scale histogram normalization method:
%                      'off'   - raw dwell counts, no normalization.
%                      'state' - each state; sum of each histogram=1. (default)
%                      'file'  - all states; sum of all histograms per file=1.
%                      'time'  - dwell counts per second of observation time
%
%      'model':        QubModel object for displaying expected distributions
%
%   See also: lifetime_exp, loadDwelltimes, removeBlinks.

%   Copyright 2007-2018 Cornell University All Rights Reserved.



%% Process input arguments

% Extract first argument graphics handle (figure or axes)
hFig = [];
ax = [];
oneAx = false;      %if true, collapse all plots into one axes.
soloWindow = true;  %if true, add menu controls for stand-alone execution

if nargin>=1 && all(ishandle(varargin{1}))
    if strcmpi(get(varargin{1},'type'), 'figure')
        hFig = varargin{1};
    elseif strcmpi(get(varargin{1},'type'), 'axes')
        ax = varargin{1};
        hFig = ancestor(ax,'figure');
        oneAx = isscalar(ax);
        cla(ax,'reset');
    else
        error('Invalid first input argument');
    end
    varargin = varargin(2:end);
    soloWindow = false;
else
    % Display results in a new figure if no outputs requested.
    if nargout==0, hFig=figure; end
end

[varargout{1:nargout}] = deal([]);
set(hFig,'pointer','watch'); drawnow;


% Prompt user for file names if not given.
if numel(varargin)<1,
    dwtfilename = getFiles('*.dwt','Choose dwell-time files');
else
    dwtfilename = varargin{1};
end
if ischar(dwtfilename), dwtfilename={dwtfilename}; end
dwtfilename = findDwt(dwtfilename,'raiseError');
if numel(dwtfilename)==0, return; end


% Assign default parameter values
persistent params;

if isempty(params)
    %FIXME: these should be defined in cascadeConstants?
    params.logX = true;
    params.dx = 0.4;
    params.removeBlinks = true;
    params.normalize = 'state';
end
params.model = [];

% Merge options, giving the user's options precedence.
if numel(varargin)>1,
    params = mergestruct( params, varargin{2} );
end

% Check parameters
if ~all(ismember(params.normalize,{'off','none','state','file','time'})),
    error('Invalid normalization option ''%s''',params.normalize);
end



%% Get dwell times from file

nFiles = numel(dwtfilename);
dwells  = cell(nFiles,1);  %consolidated list of dwell times in each state
sampling = zeros(nFiles,1);

for i=1:nFiles,
    % Load dwell-times a list per state, concatinating dwells from all traces,
    % ignoring the zero-FRET state if applicable.
    if params.removeBlinks,
        [dwells{i},sampling(i)] = loadDwelltimes( dwtfilename{i}, 'removeBlinks' );
    else
        [dwells{i},sampling(i)] = loadDwelltimes( dwtfilename{i} );
    end
end

% Verify all input idealizations are roughly consistent.
if ~all(sampling==sampling(1)),
    warning('Mismatched time resolution of input files');
end
sampling = sampling(1)/1000;  %convert to seconds.

nClasses = cellfun(@numel, dwells);
if ~all(nClasses==nClasses(1)),
    warning('Idealization models have a different number of states');
end
nClasses = max(nClasses);


% Get dwell time limits for setting axes limits later.
maxTime = 0;  %longest dwell in seconds
totalTime = zeros(nFiles,nClasses);
meanTime  = zeros(nFiles,nClasses);  %mean dwell-time per state/file.

for i=1:nFiles,
    dwellc = dwells{i};
    maxTime = max( maxTime, max(vertcat(dwellc{:})) );
    totalTime(i,:) = cellfun(@sum, dwellc)';
    meanTime(i,:)  = cellfun(@mean, dwellc)';
end



%% Calculate dwell time bins (EDGES)
if ~params.logX,
    % Linear X-axis in seconds.
    dwellaxis = 0:sampling:maxTime;
    fitaxis = dwellaxis;
else
    % Create a log time axis with a fixed number of bins.
    % histcounts uses bin edges of E(k) <= X(i) < E(k+1).
    dwellaxis = log(sampling):params.dx:log(maxTime*3);
    fitaxis = (-5:0.1:5)';  %fine-grained axis for theoretical fit curves.
    
    % Force the bins edges to be exact intervals of the time resolution.
    % The histogram will better sample the discrete nature of the data.
    maxFrames = ceil(maxTime*3/sampling);
    fullaxis = log( (1:maxFrames)*sampling )';
    dwellaxis = unique( nearestBin(dwellaxis, fullaxis) );
%     dwellaxis = unique( floor(fullaxis/sampling)*sampling );
    
    % Normalization factor to account for varying-sized bins.
    dlx = dwellaxis(2:end) - dwellaxis(1:end-1);
    dlx = [dlx dlx(end)];
end



%% Calculate histograms and optional fit lines
histograms = cell(nFiles,nClasses);

for file=1:nFiles,
    ndwells = cellfun(@numel,dwells{file});
            
    for state=1:nClasses,
        % Small constant ensures dwells fall in the correct histogram bin.
        dwellc = dwells{file}{state} +sampling/10;
        
        % Make linear-scale survival plot.
        if ~params.logX,
            counts = histc( dwellc, dwellaxis );
            histdata = sum(counts) - cumsum(counts);
            histdata = histdata/histdata(1);
        
        % Make log-scale Sine-Sigworth plot (linear ordinate).
        else
            counts = histc( log(dwellc)', dwellaxis );
            histdata = counts./dlx;  %normalize by log-space bin size
            histdata = histdata/sum(histdata);  %normalize to 1
            
            switch params.normalize
                case {'none','off'}  %raw dwell counts
                    normFact = ndwells(state);
                case 'state'  %fraction of counts in each bin for this state
                    normFact = 100;
                case 'file'  %fraction of counts in each bin across entire file
                    normFact = 100*ndwells(state)/sum(ndwells);
                case 'time'  %fraction of dwells per total observation time
                    normFact = ndwells(state)/sum(totalTime(file,:));
            end
            histdata = normFact * histdata;
        end
        
        histograms{file,state} = to_col(histdata);
    end
end



%% Calculate distribution fit lines
% Parameter values are provided by batchKinetics (final optimized model).
% FIXME: no corrections are made for missed events!
fits = zeros( numel(fitaxis), nClasses );

if isfield(params,'model') && ~isempty(params.model)
    
    m  = params.model;
    Q  = m.calcQ();               %Normalized rate matrix (rows sum to zero)
    p0 = m.calcEquilibriumP0();   %State probabities at equilibrium
    
    for class=1:nClasses,
        inclass = m.class==class;  %true for states in the same class as the current state
        Na = sum(inclass);         %number of states in this class

        % Time constants are the eigenvalues of the submatrix of Q that
        % contains all states in the current class (Qaa). See pg. 603.
        Qaa = Q(inclass,inclass);
        [X,lambda] = eig(-Qaa);
        tau = 1./diag(lambda);
        
        % Pre-exponential terms (sum to 1)
        Qfa = Q(~inclass,inclass); %rates into current class
        phi0 = ( p0(~inclass) * Qfa )  ./  ( p0(~inclass) * Qfa * ones(Na,1) );
        Y = X^-1;
        a = zeros(Na,1);
        for i=1:Na
            A = X(:,i) * Y(i,:);  %spectral matrices. See eq. 89, pg. 616
            a(i) = -tau(i) * phi0 * A * Qaa * ones(Na,1);  %eq. 41, pg. 604
        end

        % Calculate exponential mixture distribution
        expdist = 0;
        for i=1:Na
            if params.logX
                % Log scale: Sigowrth & Sine (1987) Biophys. J 52, p. 1047.
                z = fitaxis - log( tau(i) );
                expdist = expdist + a(i) * exp( z - exp(z) );
            else
                % Linear exponential decay
                expdist = expdist + a(i) * exp( -fitaxis ./ tau(i) );
            end
        end
        
        % Scale to data. FIXME
        fits(:,class) = max(histograms{1,class}) * expdist/max(expdist);
    end
end

if params.logX,
    dwellaxis = exp(dwellaxis);
    fitaxis = exp(fitaxis);
end

% Remove dark state
if params.removeBlinks
    histograms(:,1) = [];
    fits(:,1) = [];
end


% If output requested, just return the histogram data.
% If no output requested, display the histograms instead.
% If an axis is given and output is requested, do both. FIXME
if nargout>0
    output = {dwellaxis,histograms,fits};
    [varargout{1:nargout}] = output{1:nargout};
    return;
end



%% Display plots, one state per panel.

% Find a good zoom axis range for viewing all of the histograms.
h = [histograms{:}];
ymax = 1.1*max(h(:));

if params.logX,
    xmax = dwellaxis(end);
else
    xmax = dwellaxis(  find( sum(h>0.01,2), 1, 'last' )  );
end

% Choose ordinate label based on normalization
if params.logX
    switch params.normalize
        case {'none','off'}
            ordinate = 'Counts';
        case 'state'
            ordinate = 'Counts (%)';
        case 'file'
            ordinate = 'Counts (% of file)';
        case 'time'
            ordinate = 'Counts s^{-1}';
        otherwise
            error('Invalid normalization setting');
    end
else
    ordinate = 'Dwell Survival (%)';
end


% If expected mean dwell times provided, show them as fit lines.
% Here, calculate normalization constants and change histogram line style.
if isfield(params,'model') && ~isempty(params.model)
    lineStyle = '.';
else
    lineStyle = '-';
end

% Use QuB model colors for state lines for consistency.
colors = QubModelViewer.colors;
if params.removeBlinks, colors(1)=[]; end

% Draw dwell-time histograms for each state
nStates = size(histograms,2);
for state=1:nStates
    
    % Establish axes to plot in for all possible input choices
    if oneAx
        if isempty(ax)
            ax = axes('Parent',hFig);
        end
        curAx = ax;
    else
        if numel(ax) < state
            ax(state) = subplot( nStates, 1, state, 'Parent',hFig ); %#ok<AGROW>
        end
        curAx = ax(state);
    end

    % Plot dwell times as Sine-Sigworth (log scale) or surfival plot (linear)
    if params.logX,
        semilogx( curAx, dwellaxis, [histograms{:,state}], lineStyle, 'Color',colors{state} );
        set( curAx,'XTick',10.^(-4:4) )
    else
        plot( curAx, dwellaxis, [histograms{:,state}], lineStyle, 'Color',colors{state} );
    end
    hold( curAx, 'on' );

    % Draw fit lines, if applicable. (skip legend by setting HandleVisibility)
    if isfield(params,'model') && ~isempty(params.model)
        plot( curAx, fitaxis, fits(:,state), '-', 'Color',colors{state}, 'HandleVisibility','off' );
    end
    
    ylabel(curAx, ordinate);
    xlim(  curAx, [dwellaxis(1) xmax] );
    ylim(  curAx, [0 ymax] );
    if ~oneAx
        title( curAx, sprintf('State %d',state) );
        hold( curAx, 'off' );
    end
end


if oneAx
    lines = strsplit(sprintf('State #%d_',1:nStates),'_');
    legend( ax, lines(1:end-1) );
else
    legend( ax(end), trimtitles(dwtfilename) );
    linkaxes(ax,'xy');
end
xlabel( ax(end), 'Dwell Time (s)' );
set(hFig,'pointer','arrow'); drawnow;

if ~soloWindow, return; end



%% Add menu items for adjusting settings and saving output to file
prompt = {'Remove blinks:', 'Log scale:', 'Log bin size:', 'Normalization:'};
fields = {'removeBlinks', 'logX', 'dx', 'normalize'};
types{4} = {'none','state','file','time'};
cb = @(~,~)settingdlg(params,fields,prompt,types,@dwellhist,{hFig,dwtfilename});
output = [to_col(dwellaxis) horzcat(histograms{:})];

defaultFigLayout( hFig, @(~,~)dwellhist(getFiles('*.dwt'),params), ...
                      @(~,~)dwellhist(hFig,getFiles('*.dwt'),params), ...
                      {@exportTxt,dwtfilename,output}, ...
       {'Change settings...',cb; ...
        %'Reset settings',@(~,~)dwellhist(hFig,dwtfilename) ...  %FIXME!
        'Copy output',{@clipboardmat,output}}  );


    
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
function exportTxt(~,~,files,output)
% Callback function for the "save histograms" button in the histogram figure.

names = trimtitles(files);
nFiles = numel(names);
nStates = (size(output,2)-1)/nFiles;

% Ask the user for an output filename.
[f,p] = uiputfile('*.txt','Select output filename',[mfilename '.txt']);
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
dlmwrite(outFilename, output, 'delimiter','\t', '-append');

end



