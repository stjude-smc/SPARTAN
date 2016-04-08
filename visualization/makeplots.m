function varargout = makeplots(dataFilenames, titles, varargin)
%MAKEPLOTS  Creates an array of contour, pop hist, and TD plots
% 
%   MAKEPLOTS(FILES, TITLES, OPTIONS)
%   Creates a 3xN panel of 2D and 1D population histograms and TD plots,
%   where FILES specifies the locations of the datasets to plot.  TITLES
%   specifies the titles to plot above each dataset (optional).
%   If no FILES are specified, user will be asked to select them.
%   
%   OPTIONS is a structure with settings for how to display the data and
%   how to calculate histograms, etc. See cascadeConstants.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


narginchk(0,3);
nargoutchk(0,2);

updateSpartan; %check for updates

% Get default or last makeplots settings.
defaults = mpdefault();


%% Process input data and parameters

% If no files are given, prompt the user.
if nargin<1,
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
if nargin<2,
    titles = trimtitles(f);
elseif ischar(titles),
    titles = {titles};
end


% Other options (passed as a structure)
if nargin>=3,
    defaults = catstruct( defaults, varargin{1}, 'sorted' );
end

defaults.contour_bounds = [1 defaults.contour_length defaults.fretRange];



%% Display all plots

if ~isfield(defaults,'targetAxes')
    c = cascadeConstants();
    h1 = figure('Name', [mfilename ' - ' c.software]);
else
    h1 = get(defaults.targetAxes{1,1},'parent');
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
switch nargout,
    case 1
        varargout = {h1};
    case 2
        varargout = {h1,dataFilenames};
end


% Give a warning if using some funky normalization.
if isfield(defaults,'cplot_normalize_to_max') && defaults.cplot_normalize_to_max
    disp('NOTE: these plots are normalized to the plot with the largest number of traces!!');
end


% =============== ADD GUI CONTROLS ================
hMenu = findall(h1,'tag','figMenuUpdateFileNew');
delete(allchild(hMenu));
set(hMenu, 'Callback', @(~,~)makeplots() );

% hMenu = findall(h1,'tag','figMenuOpen');
% set(hMenu, 'Callback', @(~,~)makeplots(h1,getFiles(),titles,defaults) );

hMenu = findall(h1,'tag','figMenuGenerateCode');
set(hMenu, 'Label','Export .txt files', 'Callback',@(x,y)saveFiles2(x,y));
       
hEditMenu = findall(h1, 'tag','figMenuEdit');
delete(allchild(hEditMenu));
uimenu('Label','Change settings...', 'Parent',hEditMenu, 'Callback',@(x,y)changeDisplaySettings2(x,y));
uimenu('Label','Reset settings', 'Parent',hEditMenu, 'Callback',@(x,y)resetSettings2(x,y));



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
opt = settingsDialog(handles.options,fields,prompt);

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

if isempty(mpd),
    const = cascadeConstants();
    mpd = const.defaultMakeplotsOptions;
    mpd.contour_bounds = [1 mpd.contour_length mpd.fretRange];
end

if isfield(mpd,'targetAxes'),
    mpd = rmfield(mpd,'targetAxes');  %left over from rtdTool
end

output = mpd;

end %FUNCTION mpdefault




%% ===================== LOOP OVER EACH DATA FILE ====================== 

function plotData(hObject,handles)
% This is the function that actually loads the data, calculates the plots,
% and displays them.

set(handles.hFig,'pointer','watch'); drawnow;

options = handles.options;
dataFilenames = handles.dataFilenames;
baseFilenames = handles.baseFilenames;
nFiles = numel(baseFilenames);
titles = handles.titles;

N = cellfun(@sizeTraces, dataFilenames);

% Determine which files have dwell-time data so that TD plots should be
% made. This determines the number of subplot rows ahead of time.
dwtfnames = strcat(handles.baseFilenames,'.qub.dwt');
has_dwt = true(nFiles,1);

for i=1:nFiles,
    if ~exist(dwtfnames{i},'file'),
        dwtfnames{i} = [baseFilenames{i} '.dwt'];
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

[histmax,cplotax,histax,tdax] = deal( zeros(nFiles,1) );
[cplotdataAll,cpdataAll,shistAll,tdpAll,idl] = deal( cell(nFiles,1) );

% Load cached idealization data from previous calls if possible.
% FIXME: not implemented because of memory requirements.
% if isfield(handles,'idl'),
%     idl = handles.idl;
% end


% Remove latent annotations, if any.
delete(findall(handles.hFig,'Tag','Nmol'))



%%
for k=1:nFiles,
    
    % Load FRET data
    if ~exist(dataFilenames{k}, 'file')
        disp('FRET Data missing, skipping.'); continue;
    end
    
    data = loadTraces( dataFilenames{k} );
    
    % The data may have multiple FRET signals to analyze. Ask which to use.
    if data.isChannel('fret2'),
        a = questdlg('This data has multiple FRET channels. Which should be used?', ...
                     'Select FRET channel to use','fret','fret2','Cancel', 'fret');
        if strcmp(a,'Cancel'), return; end
        fret = data.(a);
        options.fretField = a;
    else
        fret = data.fret;
    end
    
    if isempty(fret),
        disp('FRET Data missing, skipping.'); continue;
    end
    
    
    %% ============== DRAW POPULATION CONTOUR HISTOGRAMS ============== 
    
    %---- LOAD OR GENERATE FRET CONTOUR PLOT DATA
    cplotdata = makecplot( fret, options );

    if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
        cplotdata(2:end,2:end) = cplotdata(2:end,2:end).*(N(k)/max(N));
    end
    cplotdataAll{k} = cplotdata;
    
    
    %---- DRAW FRET CONTOUR PLOT ----
    if ~isfield(options,'targetAxes')
        ax = subplot( nrows, nFiles, k, 'Parent',handles.hFig );
    else
        ax = options.targetAxes{k,1};
    end
    cplotax(k) = ax;
    
    % Draw the contour plot (which may be time-binned)
    % FIXME: we are passing the default options
    cpdataAll{k} = cplot( ax, cplotdata, options.contour_bounds, options );
    
    % Formatting
    title( titles{k}, 'Parent',ax );
    
    if ~options.hideText,
        ap = get(ax,'Position');  %left bottom width height
        annotation( handles.hFig, 'textbox', [ap(1)+0.925*ap(3) ap(2)+0.925*ap(4) 0.1*ap(3) 0.1*ap(4)], ...
                    'String',sprintf('N=%d', N(k)), 'HorizontalAlignment','right', ...
                    'LineStyle','none', 'tag','Nmol' );
    end
    
    if k>1,
        ylabel(ax,'');
        xlabel(ax,'');
    end
    
    
    
    %% ================ DRAW STATE OCCUPANCY HISTOGRAMS ================ 
   
    if has_dwt(k) && isempty(idl{k}),
        [dwt,~,offsets] = loadDWT(dwtfnames{k});
        idl{k} = dwtToIdl(dwt, offsets, data.nFrames, data.nTraces);
    end
    
    %---- GENERATE STATE OCCUPANCY HISTOGRAMS
    if ~options.no_statehist && has_dwt(k),
        shist = statehist( idl{k}, data, options );        
        shistAll{k} = shist;
    end
    
    
    if ~isfield(options,'targetAxes')
        ax = subplot( nrows, nFiles, nFiles+k, 'Parent',handles.hFig );
    elseif size(options.targetAxes,2)>1
        ax = options.targetAxes{k,2};
    else
        continue;
    end
    histax(k) = ax;
    cla(ax);  hold(ax,'on');
    set(ax,'ColorOrder',options.colors);
    
    
    %---- GENERATE STATE OCCUPANCY HISTOGRAMS
    % ...if no dwell-time info is available.
    if ~has_dwt(k) || options.no_statehist || isempty(shist)
        fretaxis = cplotdata(2:end,1);      
        histdata = cplotdata(2:end,2:options.contour_length+1)*100;
        pophist = nansum(histdata,2)/options.contour_length;   %normalization
        
        histmax(k) = max(max( pophist(fretaxis>0.05) ));
        
        bar( ax, fretaxis, pophist );
        

    %---- GENERATE STATE OCCUPANCY HISTOGRAMS
    % ...if dwell-time info is available.
    else
        bins = shist(:,1);
        histdata = shist(:,2:end)*100;
        [~,nStates] = size(histdata);
        
        % Pad with empty bins for display
        df = mean(diff(bins));
        bins = [bins(1)-df; bins; bins(end)+df];  %#ok
        histdata = [zeros(1,nStates); histdata; zeros(1,nStates)];  %#ok
        
        % If requested, remove 0-FRET peak and renormalize
        if options.ignoreState0
            nStates = nStates-1;
            
            histdata = histdata(:,2:end);
            histdata = 100*histdata ./ sum(histdata(:));
        end
        
        % If the option is set, rescale so that plots with only a few
        % molecules show low occupancy in the statehist.
        if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
            histdata = histdata.*(N(k)/max(N));
        end

        % Draw translucent, filled area underneath curves
        for j=1:nStates
            patch( bins, histdata(:,j), options.colors(j,:), ...
                    'EdgeColor','none','FaceAlpha',0.25, 'Parent',ax );
        end
        
        % Draw occupancy histograms as solid lines
        plot( ax, bins, histdata, 'LineWidth',1.5 );
        
        % Add a line with total occupancy.
        totalHist = sum( histdata, 2 );
        plot( ax, bins, totalHist, 'k-', 'LineWidth',1.5 );
        
        histmax(k) = max( totalHist(bins>0.05) );
    end
    
    drawnow;
    
    
    
    %% ========================= DRAW TD PLOTS ========================= 
    
    if ~has_dwt(k) || options.no_tdp
        continue;
    end
    
    
    %---- LOAD TDPLOT DATA ----
    tdp = tdplot(idl{k}, data, options);
    if tdp(1,1)==0,
        disp('Skipping invalid TDP');
        continue;
    end
    
    % If total time (normalization factor) is saved in TD plot, use it to
    % get raw number of transitions
    total_time = tdp(1,1);
    t = sum(sum( tdp(2:end,2:end)*total_time ));
    
    % If the option is set, rescale so that plots with only a few
    % molecules show low occupancy in the statehist.
    if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
        tdp(2:end,2:end) = tdp(2:end,2:end) * N(k)/max(N);
    end
    
    tdpAll{k} = tdp;
    
    
    %---- DISPLAY TD PLOT ----
    if ~isfield(options,'targetAxes')
        tdax(k) = subplot( nrows, nFiles, 2*nFiles+k, 'Parent',handles.hFig );  
    elseif size(options.targetAxes,2)>2
        tdax(k) = options.targetAxes{k,3};
    else
        continue;
    end
    tplot( tdax(k), tdp, options );
    
    % Formatting
    if ~options.hideText,
        textOpt = {'HorizontalAlignment','center', 'Parent',tdax(k)};
        text( 0.45,0.9, sprintf('N_t=%.0f',t),            textOpt{:} );
        text( 0.45,0.0, sprintf('t/s=%.2f',t/total_time), textOpt{:} );
    end
    
    drawnow;
    
end  %for each file.

set(handles.hFig,'pointer','arrow'); drawnow;



%% =================== FINISH UP =================== 

% Save results back to handles object
handles.cplotdataAll = cplotdataAll;
handles.cpdataAll = cpdataAll;
handles.shistAll = shistAll;
handles.tdpAll = tdpAll;
% handles.idl = idl;

guidata(hObject,handles);


% Finish formatting all plots, including linking the axes
if any(cplotax)
    cplotax = cplotax(cplotax~=0);
    linkaxes( cplotax, 'xy' );
    ylim( cplotax(1), options.fretRange );
    
    ylabel(cplotax(1),'FRET');
    xlabel(cplotax(1),'Time (frames)');
    set(cplotax(2:end),'YTickLabel',[]);
    
    set(cplotax, 'YGrid','on', 'Box','on');
end

if any(histax),
    histax = histax(histax~=0);
    linkaxes( histax, 'xy' );
    ylim( histax(1), [0 max(histmax)+0.5] );
    xlim( histax(1), options.fretRange );
    
    ylabel( histax(1),'Occupancy (%)' );
    xlabel( histax(1),'FRET' );
    set(histax(2:end),'yticklabel',[]);
    
    set(histax, 'YGrid','on', 'Box','on');
end

if any(tdax),
    tdax = tdax(tdax~=0);
    linkaxes(tdax, 'xy');
    ylim(tdax(1), options.fretRange );
    xlim(tdax(1), options.fretRange );
    
    ylabel(tdax(1),'Final FRET');
    xlabel(tdax(1),'Inital FRET');
    set(tdax(2:end),'yticklabel',[]);
    
    set(tdax, 'XGrid','on', 'YGrid','on', 'Box','on');
end



end  % function plotData








