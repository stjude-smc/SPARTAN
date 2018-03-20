function varargout = makeplots(varargin)
%MAKEPLOTS  Creates an array of contour, pop hist, and TD plots
% 
%   MAKEPLOTS(FILES, TITLES, OPTIONS)
%   Creates a 3xN panel of 2D and 1D population histograms and TD plots,
%   where FILES specifies the locations of the datasets to plot.  TITLES
%   specifies the titles to plot above each dataset (optional).
%   If no FILES are specified, user will be asked to select them.
%   
%   MAKEPLOTS(H,...) draws the plot in the figure H.
%   
%   OPTIONS is a structure with settings for how to display the data and
%   how to calculate histograms, etc. See cascadeConstants.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


narginchk(0,4);
nargoutchk(0,2);

updateSpartan; %check for updates

% Get default or last makeplots settings.
defaults = mpdefault();


%% Process input data and parameters

% First argument may be a target figure handle to draw in.
if nargin>=1 && isgraphics(varargin{1}(1)),
    if strcmpi(get(varargin{1},'Type'),'figure'),
        h1 = varargin{1};
        varargin(1) = [];
    else
        error('Invalid figure target');
    end
else
    h1 = [];
end
% dataFilenames, titles, varargin

% If no files are given, prompt the user.
if numel(varargin)>=1,
    dataFilenames = varargin{1};
else
    dataFilenames = getFiles('*.traces','Choose a traces file:');
end

if ~iscell(dataFilenames),
    dataFilenames = {dataFilenames};
end

nFiles = numel(dataFilenames);
if nFiles==0, return; end


% Generate plot titles and base filenames if none given.
[p,f] = cellfun(@fileparts, dataFilenames, 'UniformOutput',false);
baseFilenames = fullfile(p,f);
if numel(varargin)>=2,
    titles = varargin{2};
    if ischar(titles),
        titles = {titles};
    end
else
    titles = trimtitles(dataFilenames);
end


% Other options (passed as a structure)
if numel(varargin)>=3,
    defaults = orderfields(  mergestruct( defaults, varargin{3} )  );
end

defaults.contour_bounds = [1 defaults.contour_length defaults.fretRange];



%% Display all plots

if isempty(h1),
    if ~isfield(defaults,'targetAxes')
        c = cascadeConstants();
        h1 = figure('Name', [mfilename ' - ' c.software]);
    else
        h1 = get(defaults.targetAxes{1,1},'parent');
    end
else
    clf(h1);
end

% Force the cursor to be normal arrow. When loading from .fig files, it is
% sometimes mysteriously 'watch'. MATLAB bug?
set(h1,'CreateFcn','set(gcbo,''pointer'',''arrow'')');

newHandles.options = defaults;
newHandles.baseFilenames = baseFilenames;
newHandles.dataFilenames = dataFilenames;
newHandles.titles = titles;
newHandles.hFig = h1;

plotData(h1,newHandles);


% Set outputs, if requested.
output = {h1,dataFilenames};
[varargout{1:nargout}] = output{1:nargout};


% =============== ADD GUI CONTROLS ================
defaultFigLayout( h1, @(~,~)makeplots, @(~,~)makeplots(h1), @(x,y)saveFiles2(x,y), ...
                  {'Change settings...',@(x,y)changeDisplaySettings2(x,y); ...
                   'Reset settings',    @(x,y)resetSettings2(x,y)} );


% =============== COMPATIBILITY with 2.11 and earlier ================
% Display a warning instead of crashing.
function saveFiles(varargin) %#ok<DEFNU>
    disp('This feature is not available for versions before 2.12.');
end

function changeDisplaySettings(varargin) %#ok<DEFNU>
    disp('This feature is not available for versions before 2.12.');
end

function resetSettings(varargin) %#ok<DEFNU>
    disp('This feature is not available for versions before 2.12.');
end



end %FUNCTION MAKEPLOTS




%% ================ GUI CALLBACKS ================ 
% Called when any of the buttons at the bottom of the makeplots window are
% clicked. All plot data is stored in guidata.

function changeDisplaySettings2(hObject,~)
% Prompt user to change any makeplots display settings.

handles = guidata(hObject);

% 1. Get the new value from the user.
prompt = {'Contour length (frames):', 'Contour offset (frames):', ...
          'Contour scaling factor:', 'FRET bin size:', ...
          'Hide photobleaching:', 'TD plot scaling factor:', ...
          'Hide blinks in TD plots', 'Truncate TD plot'};
      
fields = {'contour_length', 'pophist_offset', ...
          'cplot_scale_factor', 'contour_bin_size', ...
          'cplot_remove_bleached', 'tdp_max', ...
          'hideBlinksInTDPlots', 'truncate_tdplot' };
opt = settingdlg(handles.options,fields,prompt);
if isempty(opt), return; end

opt.fret_axis = -0.1:opt.contour_bin_size:1.2;
opt.contour_bounds = [1 opt.contour_length opt.fretRange];
handles.options = opt;
mpdefault(opt);  %save as starting values for future calls to makeplots.

% 3. Redraw plots.
plotData(hObject,handles);

end %FUNCTION changeDisplaySettings



function resetSettings2(hObject,~)
% Reset all display settings to their defaults in cascadeConstants.

handles = guidata(hObject);
handles.options = mpdefault([]);  %reset to default, save persistent state.
plotData(hObject,handles);

end %FUNCTION resetSettings



function saveFiles2(hObject,~)
% Save plot data to text files for importing and plotting in Origin.
    
handles = guidata(hObject);
base = handles.baseFilenames;

for k=1:numel(base),
    % Save contour FRET histogram
    if ~isempty(handles.cplotdataAll{k}),
        dlmwrite( [base{k} '_hist.txt'] , handles.cplotdataAll{k}, ' ' );
    end

    % Save display-format (time-binned) contour FRET histogram.
    if ~isempty(handles.cpdataAll{k}),
        dlmwrite( [base{k} '_displayhist.txt'], handles.cpdataAll{k}, ' ');
    end

    % State histogram
    if ~isempty(handles.shistAll{k}),
        dlmwrite( [base{k} '.qub_shist.txt'], handles.shistAll{k}, ' ' );
    end

    % TD Plot
    if ~isempty(handles.tdpAll{k}),
        dlmwrite( [base{k} '.qub_tdp.txt'], handles.tdpAll{k}, ' ' );
    end
end
    
end %FUNCTION saveFiles


function output = mpdefault(input)
% Read and write access to persistent makeplots settings.
% mpdefault(INPUT) saves the INPUT state. empty input resets to defaults.
% OUTPUT = mpdefault() get the current persistent state.

persistent mpd;

if nargin>0, mpd=input; end

% Reload defaults if persistent variable is not set, or
% cascadeCosntants.m has been modified (timestamp has changed).
const = cascadeConstants();
if isempty(mpd) || mpd.tstamp~=const.tstamp,  
    mpd = const.defaultMakeplotsOptions;
    mpd.tstamp = const.tstamp;
    mpd.contour_bounds = [1 mpd.contour_length mpd.fretRange];
end

if isfield(mpd,'targetAxes'),
    mpd = rmfield(mpd,'targetAxes');  %left over from rtdTool
end

output = mpd;

end %FUNCTION mpdefault




%% ===================== LOOP OVER EACH DATA FILE ====================== 

function plotData(hObject,handles)
% Load data, calculate plots, and display them.
set(handles.hFig,'pointer','watch'); drawnow;

options = handles.options;
dataFilenames = handles.dataFilenames;
nFiles = numel(dataFilenames);

N = cellfun(@sizeTraces, dataFilenames);
if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
    options.cplot_normalize_to_max = max(N);
    disp('NOTE: plots are normalized to largest number of traces!');
end

% Determine which files have dwell-time data so that TD plots should be
% made. This determines the number of subplot rows ahead of time.
dwtfnames = strcat(handles.baseFilenames,'.qub.dwt');
has_dwt = true(nFiles,1);

for i=1:nFiles,
    if ~exist(dwtfnames{i},'file'),
        dwtfnames{i} = [handles.baseFilenames{i} '.dwt'];
    end
    
    if ~exist(dwtfnames{i},'file'),
        dwtfnames{i} = [];
        has_dwt(i) = false;
    end
end

if any(has_dwt) && ~options.no_tdp, 
    nrows = 3;
else
    nrows = 2;
end

[idl,histmax,cplotax,histax,tdax] = deal( zeros(nFiles,1) );
[cplotdataAll,cpdataAll,shistAll,tdpAll] = deal( cell(nFiles,1) );


% Remove latent annotations, if any.
delete(findall(handles.hFig,'Tag','Nmol'))


% Get target axes for plots.
% FIXME this won't handle empty targetAxes elements correctly!
if isfield(options,'targetAxes')
    cplotax = [options.targetAxes{:,1}];
    
    if size(options.targetAxes,2) >= 2,
        histax  = [options.targetAxes{:,2}];
    end
    if size(options.targetAxes,2) >= 3,
        tdax    = [options.targetAxes{:,3}];
    end
    
else
    axopt = {'Parent',handles.hFig};
    for k=1:nFiles,
        cplotax(k) = subplot(nrows, nFiles, k, axopt{:});
        histax(k)  = subplot(nrows, nFiles, nFiles+k, axopt{:});
        if nrows>2 && has_dwt(k),
            tdax(k) = subplot(nrows, nFiles, 2*nFiles+k, axopt{:});
        end
    end
end


% Calculate and display all plots for each file:
options.fretFieldLabel = 'FRET';

for k=1:nFiles,
    
    % Load FRET data
    if ~exist(dataFilenames{k}, 'file')
        disp('FRET Data missing, skipping.'); continue;
    end
    
    data = loadTraces( dataFilenames{k} );
    
    % The data may have multiple FRET signals to analyze. Ask which to use.
    if k==1 && data.isChannel('fret2'),
        a = questdlg('This data has multiple FRET channels. Which should be used?', ...
                     'Select FRET channel to use','fret','fret2','Cancel', 'fret');
        if strcmp(a,'Cancel'), return; end
        options.fretField = a;
        options.fretFieldLabel = [ data.fretAxisLabel num2str(strcmpi(a,'fret2')+1) ];
    end
    
    
    %% ============== POPULATION CONTOUR HISTOGRAMS ============== 
    % Display contour plots. (cpdataAll is truncated and time-binned)
    [cplotdataAll{k},cpdataAll{k}] = cplot(cplotax(k), data, options);
    title( handles.titles{k}, 'Parent',cplotax(k) );
    
    if ~options.hideText,
        ap = get(cplotax(k),'Position');  %left bottom width height
        annotation( handles.hFig, 'textbox', [ap(1)+0.925*ap(3) ap(2)+0.925*ap(4) 0.1*ap(3) 0.1*ap(4)], ...
                    'String',sprintf('N=%d', N(k)), 'HorizontalAlignment','right', ...
                    'LineStyle','none', 'tag','Nmol' );
    end
    
    if k>1,
        ylabel(cplotax(k), '');
        xlabel(cplotax(k), '');
    end
    
    
    %% ================ STATE OCCUPANCY HISTOGRAMS ================ 
    if has_dwt(k),
        try
            [dwt,~,offsets] = loadDWT(dwtfnames{k});
            idl = dwtToIdl(dwt, offsets, data.nFrames, data.nTraces);
        catch e
            disp(['Invalid idealization: ' e.message]);
            has_dwt(k) = false;    
        end
    end
    
    if histax(k)==0, continue; end
    
    %---- FRET histogram for if data was not idealized.
    if ~has_dwt(k) || options.no_statehist
        cplotdata = cplotdataAll{k};
        fretaxis = cplotdata(2:end,1);      
        histdata = cplotdata(2:end,2:options.contour_length+1)*100;
        pophist = nansum(histdata,2)/options.contour_length;   %normalization
        
        histmax(k) = max(max( pophist(fretaxis>0.05) ));
        barh(histax(k), fretaxis, pophist);
        
    %---- Statehist: FRET histogram for each idealized state.
    else
        [shistAll{k}, histmax(k)] = statehist(histax(k), idl, data, options);
    end
    
    set(histax(k),'ytick', 0:0.2:1);
%     drawnow;
    
        
    %% ========================= TD PLOTS ========================= 
    if nrows<3 || ~has_dwt(k) || tdax(k)==0, continue; end
    
    % Calculate TD plot
    tdp = tdplot(idl, data, options);
    if tdp(1,1)==0,
        disp('Skipping invalid TDP');
        continue;
    end
    tdpAll{k} = tdp;
    
    % Show TD plot
    tplot(tdax(k), tdp, options);
    drawnow;
    
    
end  %for each file.

set(handles.hFig,'pointer','arrow'); drawnow;



%% =================== FINISH UP =================== 

% Save results back to handles object
handles.cplotdataAll = cplotdataAll;
handles.cpdataAll = cpdataAll;
handles.shistAll = shistAll;
handles.tdpAll = tdpAll;

guidata(hObject,handles);


% Finish formatting all plots, including linking the axes
if any(cplotax~=0)
    cplotax = cplotax(cplotax~=0);
    linkaxes( cplotax, 'xy' );
    ylim( cplotax(1), options.fretRange );
    
    ylabel(cplotax(1), options.fretFieldLabel);
    xlabel(cplotax(1),'Time (s)');  %fixme could be frames
    set(cplotax(2:end),'YTickLabel',[]);
    
    set(cplotax, 'YGrid','on', 'Box','on');
end

if any(histax~=0),
    histax = histax(histax~=0);
    linkaxes( histax, 'xy' );
    xlim( histax(1), [0 max(histmax)+0.5] );
    ylim( histax(1), options.fretRange );
    
    xlabel( histax(1),'Counts (%)' );
    ylabel( histax(1), options.fretFieldLabel );
    set(histax(2:end),'yticklabel',[]);
    
    set(histax, 'YGrid','on', 'Box','on');
end

if any(tdax~=0),
    tdax = tdax(tdax~=0);
    linkaxes(tdax, 'xy');
    ylim(tdax(1), options.fretRange );
    xlim(tdax(1), options.fretRange );
    
    ylabel(tdax(1),'Final FRET');
    xlabel(tdax(1),'Inital FRET');
    set(tdax(2:end),'yticklabel',[]);
    
    set(tdax, 'XGrid','on', 'YGrid','on', 'Box','on');
end

zoom(handles.hFig,'on');


end  % function plotData








