function varargout = transitionsPerSecond(varargin)
% transitionsPerSecond  Calculate average transition rates
%
%   [TPS,ERR] = transitionsPerSecond( FILES ) calculates the average number of
%   transitions per second in each .dwt file in the cell array FILES.
%   Transitions to and from the dark state (blinking) are ignored.
%   ERR are bootstrapped standard errors for each file.
%   
%   transitionsPerSecond(FILES,PARAMS) specifies optional parameters:
%
%    - removeZeroState: If true, do not consider the zero-FRET state (class 1)
%    - truncateLength:  Number of frames at the beginning to use for calculation
%    - errorbars:       'bootstrap' or 'std' (see separateDays.m)
%
%   transitionsPerSecond(...) if no outputs requested, display in a new figure.
%
%   transitionsPerSecond(AX,...) plots in the scalar axes AX.
%
%   See also: percentTime, makeplots.

%   Copyright 2007-2017 Cornell University All Rights Reserved.


% Default parameter values
persistent params;

if isempty(params)
    params.truncateLength = Inf;
    params.removeZeroState = true;
    params.errorbars = 'bootstrap';
end


%% Process input arguments
narginchk(0,3);
nargoutchk(0,2);
[varargout{1:nargout}] = deal([]);
[cax,args] = axescheck(varargin{:});

switch numel(args)
    case 0
        dwtFilenames = getFiles('*.dwt','Choose idealization files:');
    case 1
        dwtFilenames = args{1};
    case 2
        [dwtFilenames,inputParams] = args{:};
        params = mergestruct(params, inputParams);
end
if ischar(dwtFilenames),
    dwtFilenames = {dwtFilenames};
end
nFiles = numel(dwtFilenames);
if nFiles==0, return; end

% Find associated data file types, regardless of input format.
[trcfiles,dwtFilenames] = findTracesDwt( dwtFilenames );


% Make a new figure or clear user-defined target.
if ~isempty(cax),
    hFig = ancestor(cax,'figure');
    set(hFig, 'pointer','watch'); drawnow;

elseif nargout==0,
    hFig = figure;
    set(hFig, 'pointer','watch'); drawnow;
end



%% Calculate transitions/sec for each datafile using bootstrap sampling.
bootfun = @(N,time) sum(N)/sum(time);  %transitions per second

meanTPS = zeros(nFiles,1);
stdTPS = zeros(nFiles,1);

for i=1:nFiles,    
    
    switch lower(params.errorbars)
    case {'bootstrap','none'}
        % Load DWT data
        [dwells,sampling] = loadDWT( dwtFilenames{i} );
    
        % Calculate bootstrap samples to estimate standard error.
        [nEvents,totalTime] = dwellCount(dwells, sampling, params);
        meanTPS(i,:) = sum(nEvents)/sum(totalTime);
        
        if ~strcmpi(params.errorbars,'none'),
            stdTPS(i,:) = std(  bootstrp(1000, bootfun, nEvents,totalTime)  );
        end
        
    case {'stdev','std'}
        [datarep,idlrep] = separateDays( trcfiles{i}, dwtFilenames{i} );
        reps = zeros(numel(datarep),1);
        
        for r=1:numel(datarep)
            [nEvents,totalTime] = dwellCount(idlrep{r}, datarep{r}.sampling, params);
            reps(r) = sum(nEvents)/sum(totalTime);
        end
        meanTPS(i) = mean(reps);
        stdTPS(i)  = std(reps);
        
    otherwise
        error('Unknown errorbar calculation method');
    end
    
end %for each file

% Handle output parameters
outputs = {meanTPS,stdTPS};
[varargout{1:nargout}] = outputs{1:nargout};

if nargout>0 && isempty(cax),
    return;
end



%% Display the result.
cax = newplot(hFig);

nStates = size(meanTPS,2);
if ~strcmpi(params.errorbars,'none')
    errorbar( cax, repmat(1:nFiles,nStates,1)', meanTPS, stdTPS/2 );
else
    plot( cax, repmat(1:nFiles,nStates,1)', meanTPS, 'bo-' );
end

ylabel(cax, 'Transitions per second');
xlim(cax, [0.5 nFiles+0.5]);
set(cax,'XTick',1:nFiles);

labels = trimtitles(dwtFilenames);
if max(cellfun(@numel,labels))<20,
    set(cax,'XTick',1:nFiles, 'XTickLabelRotation',30);
    set(cax,'XTickLabel',labels);
else
    xlabel(cax,'File number');
end



%% Add menus to change settings, get data, open new plots, etc.
types  = { @isscalar, [], {'bootstrap','stdev','none'} };
fields = {'truncateLength','removeZeroState','errorbars'};
prompt = {'Use only first N frames','Ignore dark state','Error bar calculation'};

defaultFigLayout( hFig,  @(~,~)transitionsPerSecond(getFiles('*.dwt'),params), ...
                         @(~,~)transitionsPerSecond(cax,getFiles('*.dwt'),params), [], ...
   {'Change settings...',@(~,~)settingdlg(params,fields,prompt,types,mfilename,@transitionsPerSecond,{cax,dwtFilenames}) ; ...
    'Reset settings',    @(~,~)transitionsPerSecond(cax,dwtFilenames) ; ...
    'Copy values',{@clipboardmat,[meanTPS stdTPS]}  }    );

set(hFig,'pointer','arrow'); drawnow;


end %FUNCTION transitionsPerSecond.




%%
function [nEvents,totalTime] = dwellCount( input, sampling, params )
% Number of state transitions and total time alive in each segment.

% If input is an idealization matrix, convert it to dwell times.
if ischar(input)
    dwells = loadDWT(input);
elseif ~iscell(input),
    dwells = idlToDwt(input);
else
    dwells = input;
end

nTraces = numel(dwells);
nEvents = zeros( nTraces, 1 );
totalTime = zeros( nTraces, 1 );

for i=1:nTraces,
    states = dwells{i}(:,1);
    times  = dwells{i}(:,2);

    if params.removeZeroState,
        % Remove dwells in lowest FRET state (assuming it is the dark state)
        times  = times(states>1);
        states = states(states>1);
    end

    if ~isempty(times),
        % Truncate idealization
        idl = dwtToIdl( [states times] );
        idl = idl( 1:min(params.truncateLength,numel(idl)) );
        newDwt = RLEncode(idl);
        times  = newDwt(:,2);

        % Calculate number of events and total time in this trace.
        nEvents(i)   = numel(times)-1;
        totalTime(i) = sum(times)*sampling/1000; %in seconds.
    end

end %for each trace

end %function tps




