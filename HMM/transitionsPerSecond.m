function varargout = transitionsPerSecond(varargin)
% transitionsPerSecond  Calculate average transition rates
%
%   [TPS,ERR] = transitionsPerSecond( FILES ) calculates the average number of
%   transitions per second in each .dwt file in the cell array FILES.
%   Transitions to and from the dark state (blinking) are ignored.
%   ERR are bootstrapped standard errors for each file.
%   
%   [...] = transitionsPerSecond(FILES,PARAMS) specifies optional parameters:
%     'removeZeroState': Do not consider the zero-FRET state (class 1).  (true)
%
%   transitionsPerSecond(...) if no outputs requested, display in a new figure.
%
%   transitionsPerSecond(AX,...) plots in the scalar axes AX.
%
%   See also: percentTime, makeplots.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


% Default parameter values
% params.truncateLength = Inf;
params.removeZeroState = true;



%% Process input arguments
narginchk(0,3);
nargoutchk(0,2);
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

% If .traces files are given, silently look for associated .dwt file.
dwtFilenames = findDwt(dwtFilenames);



%% Calculate transitions/sec for each datafile using bootstrap sampling.
bootfun = @(N,time) sum(N)/sum(time);  %transitions per second

meanTPS = zeros(0,0);
stdTPS = zeros(0,0);

for i=1:nFiles,    
    % Load DWT data
    [dwells,sampling] = loadDWT( dwtFilenames{i} );
    nTraces = length(dwells);

    nEvents = zeros( nTraces, 1 );
    totalTime = zeros( nTraces, 1 );

    for trace=1:nTraces,
        states = dwells{trace}(:,1);
        times  = dwells{trace}(:,2);

        if params.removeZeroState,
            % Remove dwells in lowest FRET state (assuming it is the dark state)
            times  = times(states>1);
            states = states(states>1);

            % Combine dwells that are now in the same state by converting into
            % an idealization and then back to a dwell-time sequence.
            if ~isempty(times),
                idl = dwtToIdl( [states times] );
                newDwt = RLEncode(idl);
                times  = newDwt(:,2);
            end
        end

        % Calculate number of events and total time in this trace.
        nEvents(trace)   = numel(times)-1;
        totalTime(trace) = sum(times)*sampling/1000; %in seconds.

    end %for each trace

    % Calculate bootstrap samples to estimate standard error.
    meanTPS(i,:) = bootfun(nEvents,totalTime);
    stdTPS(i,:) = std(  bootstrp(1000, bootfun, nEvents,totalTime)  );
    
end %for each file

% Handle output parameters
outputs = {meanTPS,stdTPS};
[varargout{1:nargout}] = outputs{1:nargout};

if nargout>0 && isempty(cax),
    return;
end



%% Display the result.

% Make a new figure or clear user-defined target.
if isempty(cax),
    hFig = figure;
else
    hFig = get(cax,'Parent');
end
cax = newplot(hFig);

nStates = size(meanTPS,2);
errorbar( cax, repmat(1:nFiles,nStates,1)', meanTPS, stdTPS/2 );

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
hMenu = findall(hFig,'tag','figMenuUpdateFileNew');
delete(allchild(hMenu));
set(hMenu, 'Callback', @(~,~)transitionsPerSecond(getFiles('*.dwt'),params) );

hMenu = findall(hFig,'tag','figMenuOpen');
set(hMenu, 'Callback', @(~,~)transitionsPerSecond(cax,getFiles('*.dwt'),params) );


hEditMenu = findall(hFig, 'tag','figMenuEdit');
delete(allchild(hEditMenu));

uimenu('Label','Copy values', 'Parent',hEditMenu, 'Callback',{@clipboardmat,[meanTPS stdTPS]});




end %FUNCTION transitionsPerSecond.





