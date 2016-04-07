function varargout = percentTime(varargin)
% PERCENTTIME  Stable state probabilities
%
%   [PT,SE] = percentTime(FILES) returns the overall occupancy in each FRET  
%   state (PT) and bootstrapped standard errors (SE) from .dwt files in the
%   cell array FILES. States are in columns, files in rows.
%
%   [...] = percentTime(FILES,PARAMS) specifies optional parameters:  (defaults)
%     'truncateLength': number of frames from beginning to use.       (Inf)
%     'hideZeroState':  Do not show zero-FRET state (class 1).        (true)
% 
%   percentTime(...) if no outputs are requested, data are displayed instead.
%   
%   percentTime(AX,...) draws the plot in the scalar axes AX.

%   Copyright 2007-2016 Cornell University All Rights Reserved.

% FIXME: if targetting other axes, this will alter the menus...
% Could consider, generally, to use zoom context menu instead...


% Default parameter values
params.truncateLength = Inf;
params.hideZeroState = true;


%% Process input arguments
narginchk(0,3);
nargoutchk(0,2);
[cax,args] = axescheck(varargin{:});

switch numel(args)
    case 0
        filenames = getFiles('*.dwt','Choose idealization files:');
    case 1
        filenames = args{1};
    case 2
        [filenames,inputParams] = args{:};
        params = mergestruct(params, inputParams);
end

% If .traces files are given, silently look for associated .dwt file.
filenames = findDwt(filenames);
nFiles = numel(filenames);
if nFiles==0, return; end


% Get axes target if given. If not and no outputs requested, make a new one.
if ~isempty(cax),
    hFig = get(cax,'Parent');
    set(hFig, 'pointer','watch'); drawnow;
end

if isempty(cax) && nargout==0,
    hFig = figure;
    set(hFig, 'pointer','watch'); drawnow;
end



%% Load dwell-times and calculate percent time in each state for each file.
bootfun = @(times) 100*sum(times)/sum(times(:));

for i=1:nFiles,
    
    % Load dwell-time information and convert to state assignment matrix.
    [dwt,~,~,model] = loadDWT(filenames{i});
    nStates = size(model,1);
    
    idl = dwtToIdl(dwt);
    [nTraces,len] = size(idl);

    % Truncate the idealization if necessary
    if nargin>=2,
        idl = idl( :, 1:min(len,params.truncateLength) );
    end

    % Calculate percent time of each trace seperately
    tracePT = zeros(nTraces, nStates);

    for state=1:nStates,
        tracePT(:,state) = sum(idl==state,2);
    end

    % Remove zero state from consideration:
    if params.hideZeroState,
        tracePT = tracePT(:,2:end);
    end

    % Calculate bootstrap samples to estimate standard error.
    if i==1,
        meanPT = zeros(nFiles, size(tracePT,2));
        stdPT  = zeros(nFiles, size(tracePT,2));
    end
    
    meanPT(i,:) = bootfun(tracePT);
    stdPT(i,:)  = std(  bootstrp(1000, bootfun, tracePT)  );
end

% Set output values
output = {meanPT,stdPT};
[varargout{1:nargout}] = output{1:nargout};
if nargout>0, return; end



%% Plot the results
nStates = size(meanPT,2);
cax = newplot(hFig);
errorbar( cax, repmat(1:nFiles,nStates,1)', meanPT, stdPT/2 );

% Construct titles with the state number and FRET values.
states = (1:nStates)+params.hideZeroState;
fret = model(states,1);

titles = cell(nStates,1);
for i=1:nStates,
    titles{i} = sprintf('State %d (%.2f)\t',states(i),fret(i));
end
legend(cax,titles);
ylabel(cax,'Fraction occupancy');

xlim(cax,[0.5 nFiles+0.5]);
set(cax,'XTick',1:nFiles);

labels = trimtitles(filenames);
if max(cellfun(@numel,labels))<20,
    set(cax,'XTick',1:nFiles, 'XTickLabelRotation',30);
    set(cax,'XTickLabel',labels);
else
    xlabel(cax,'File number');
end

set(hFig, 'pointer','arrow');  drawnow;



%% Add menus to change settings, get data, open new plots, etc.
hMenu = findall(hFig,'tag','figMenuUpdateFileNew');
delete(allchild(hMenu));
set(hMenu, 'Callback', @(~,~)percentTime(getFiles('*.dwt'),params) );

hMenu = findall(hFig,'tag','figMenuOpen');
set(hMenu, 'Callback', @(~,~)percentTime(cax,getFiles('*.dwt'),params) );


hEditMenu = findall(hFig, 'tag','figMenuEdit');
delete(allchild(hEditMenu));
cb = @(~,~) settingsDialog(params,@percentTime,{cax,filenames});
uimenu('Label','Change settings...', 'Parent',hEditMenu, 'Callback',cb);

uimenu('Label','Copy values', 'Parent',hEditMenu, 'Callback',{@clipboardmat,meanPT});
uimenu('Label','Copy errors', 'Parent',hEditMenu, 'Callback',{@clipboardmat,stdPT});



end % FUNCTION percentTime
 


