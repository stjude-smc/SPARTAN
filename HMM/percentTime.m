function varargout = percentTime(varargin)
% PERCENTTIME  Stable state probabilities
%
%   [PT,SE] = percentTime(FILES) returns the overall occupancy in each FRET  
%   state (PT) and bootstrapped standard errors (SE) from .dwt files in the
%   cell array FILES. States are in columns, files in rows.
%
%   [...] = percentTime(FILES,PARAMS) specifies optional parameters:  (defaults)
%
%    - removeZeroState: If true, do not consider the zero-FRET state (class 1)
%    - truncateLength:  Number of frames at the beginning to use for calculation
%    - errorbars:       'bootstrap' or 'std' (see separateDays.m)
% 
%   percentTime(...) if no outputs are requested, data are displayed instead.
%   
%   percentTime(AX,...) draws the plot in the scalar axes AX.

%   Copyright 2007-2016 Cornell University All Rights Reserved.

% FIXME: if targetting other axes, this will alter the menus...
% Could consider, generally, to use zoom context menu instead...


% Default parameter values
persistent params;

if isempty(params)
    params.truncateLength = Inf;
    params.removeZeroState = true;
    params.errorbars = 'bootstrap';
end


%% Process idl arguments
narginchk(0,3);
nargoutchk(0,2);
[varargout{1:nargout}] = deal([]);
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
[trcfiles,filenames] = findTracesDwt( filenames );
nFiles = numel(filenames);
if nFiles==0, return; end


% Get axes target if given. If not and no outputs requested, make a new one.
if ~isempty(cax),
    hFig = ancestor(cax,'figure');
    set(hFig, 'pointer','watch'); drawnow;

elseif nargout==0,
    hFig = figure;
    cax = newplot(hFig);
    set(hFig, 'pointer','watch'); drawnow;
end



%% Load dwell-times and calculate percent time in each state for each file.
bootfun = @(times) 100*sum(times)/sum(times(:));

meanPT = zeros(0,0);
stdPT = zeros(0,0);

for i=1:nFiles,
    [dwt,~,~,model] = loadDWT(filenames{i});
    
    switch lower(params.errorbars)
    case {'bootstrap','none'}
        pt = tracePT(dwt, params);
        meanPT(i,:) = bootfun(pt);
        
        if ~strcmpi(params.errorbars,'none')
            stdPT(i,:)  = std(  bootstrp(1000, bootfun, pt)  );
        end
        
    case {'std','stdev'}
        [~,idlrep] = separateDays( trcfiles{i}, filenames{i} );
        reps = zeros(0,0);
        
        for r=1:numel(idlrep)
            reps(r,:) = bootfun(  tracePT(idlrep{r},params)  );
        end
        meanPT(i,:) = mean(reps);
        stdPT(i,:)  = std(reps);
        
    otherwise
        error('Unknown errorbar calculation method');
    end
end

% Set output values
output = {meanPT,stdPT};
[varargout{1:nargout}] = output{1:nargout};
if isempty(cax) && nargout>0, return; end



%% Plot the results
nStates = size(meanPT,2);

if ~strcmpi(params.errorbars,'none')
    errorbar( cax, repmat(1:nFiles,nStates,1)', meanPT, stdPT/2 );
else
    plot( cax, repmat(1:nFiles,nStates,1)', meanPT, 'o-' );
end

% Construct titles with the state number and FRET values.
states = (1:nStates)+params.removeZeroState;
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



%% Add menus to change settings, get data, open new plots, etc.
types  = { @isscalar, [], {'bootstrap','stdev','none'} };
fields = {'truncateLength','removeZeroState','errorbars'};
prompt = {'Use only first N frames','Ignore dark state','Error bar calculation'};

defaultFigLayout( hFig, @(~,~)percentTime(getFiles('*.dwt'),params), ...
                        @(~,~)percentTime(cax,getFiles('*.dwt'),params), [], ...
   {'Change settings...',@(~,~)settingdlg(params,fields,prompt,types,mfilename,@percentTime,{cax,filenames}) ; ...
    'Reset settings',   @(~,~)percentTime(cax,filenames)
    'Copy values',{@clipboardmat,meanPT}   ; ...
    'Copy errors',{@clipboardmat,stdPT}    }     );

set(hFig, 'pointer','arrow');  drawnow;


end % FUNCTION percentTime
 



%%
function output = tracePT(input, params)
% Calculate fraction occupancy in each trace

% Load dwell-time information and convert to state assignment matrix.
if ischar(input)
    input = loadDWT(input);
end
if iscell(input)
    input = dwtToIdl(input);
end
[nTraces,len] = size(input);
nStates = max(input);

% Truncate the idealization if necessary
input = input( :, 1:min(len,params.truncateLength) );

% Calculate percent time of each trace seperately.
output = zeros(nTraces, 0);
for state=1:nStates,
    output(:,state) = sum(input==state,2);
end

% Remove zero state from consideration:
if params.removeZeroState,
    output = output(:,2:end);
end

end







