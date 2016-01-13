function [h1,dataFilenames] = makeplots(dataFilenames, titles, varargin)
%MAKEPLOTS  Creates an array of contour, pop hist, and TD plots
% 
%   H = MAKEPLOTS(FILES, TITLES, OPTIONS)
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


% On first run, get default parameter values from cascadeConstants. In future
% runs, use the settings as adjusted from the previous instance.
persistent defaults;
c = cascadeConstants();

if isempty(defaults),
    defaults = c.defaultMakeplotsOptions;
end



%% Process input data and parameters

% If no files are given, prompt the user.
if nargin<1,
    disp('Select traces files and hit cancel when finished');
    dataFilenames = getFiles('*.traces','Choose a traces file:');
end

if ~iscell(dataFilenames),
    dataFilenames = {dataFilenames};
end

nFiles = numel(dataFilenames);
if nFiles == 0,
    return;
end


if nargin>=2 && ischar(titles)
    titles = {titles};
end

% Generate plot titles and base filenames if none given.
baseFilenames = cell(nFiles,1);

for i=1:nFiles,
    [p,f] = fileparts( dataFilenames{i} );
    baseFilenames{i} = fullfile(p,f);

    if nargin<2 || isempty(titles),
        titles{i} = strrep(f,'_',' ');
    end
end


% Other options (passed as a structure)
if nargin>=3,
    defaults = catstruct( defaults, varargin{1}, 'sorted' );
end

defaults.contour_bounds = [1 defaults.contour_length defaults.fretRange];



%% Display all plots

if ~isfield(defaults,'targetAxes')
    h1 = figure('Name', [mfilename ' - ' c.software]);
else
    h1 = get(defaults.targetAxes{1,1},'parent');
end

newHandles.options = defaults;
newHandles.baseFilenames = baseFilenames;
newHandles.dataFilenames = dataFilenames;
newHandles.titles = titles;
newHandles.hFig = h1;

plotData(h1,newHandles);


% Give a warning if using some funky normalization.
if isfield(defaults,'cplot_normalize_to_max') && defaults.cplot_normalize_to_max
    disp('NOTE: these plots are normalized to the plot with the largest number of traces!!');
end


% =============== ADD GUI CONTROLS ================ 
% These are buttons and things to save the plots or manipulate them without
% calling makeplots again.
% Position from LL corner is defined as [left bottom width height].
% FIXME: change to normalized units so the buttons are scaled?
% NOTE: must use anonymous function reference or saved .fig files may not find 
% these functions if makeplots.m is changed.

uicontrol( 'Style','pushbutton', 'String','Save files', ...
           'Position',[20 20 80 30], 'Callback',@(x,y)saveFiles2(x,y), ...
           'Parent',h1 );

uicontrol( 'Style','pushbutton', 'String','Change settings', ...
           'Position',[110 20 130 30], 'Callback',@(x,y)changeDisplaySettings2(x,y), ...
           'Parent',h1 );
       
uicontrol( 'Style','pushbutton', 'String','Reset settings', ...
           'Position',[250 20 120 30], 'Callback',@(x,y)resetSettings2(x,y), ...
           'Parent',h1 );


% =============== COMPATIBILITY with 2.11 and earlier ================
% The implementation of these functions changed in 2.12. Before, all data were
% in the local stack so the buttons would not function correctly when loading
% the figures from file. Disabled to prevent errors or confusion.
% Without these, users cannot load old .fig files.
function saveFiles(varargin) %#ok<DEFNU>
    disp('This feature is not available for versions before 2.12.');
end

function changeDisplaySettings(varargin) %#ok<DEFNU>
    disp('This feature is not available for versions before 2.12.');
end

function resetSettings(varargin) %#ok<DEFNU>
    disp('This feature is not available for versions before 2.12.');
end



%% ================ GUI CALLBACKS ================ 
% Called when any of the buttons at the bottom of the makeplots window are
% clicked. All plot data is stored in guidata. These are nested in the main
% function scope to share persistent settings across plots.

function changeDisplaySettings2(hObject,~)
% Changes how much (how many frames) of the movie to show in plots.
% This is equivalent to changing pophist_sumlen in cascadeConstants.

handles = guidata(hObject);
opt = handles.options;

% 1. Get the new value from the user.
prompt = {'Contour length (frames):', 'Contour offset (frames):', ...
          'Contour scaling factor:', 'FRET bin size:', ...
          'Hide photobleaching:', 'TD plot scaling factor:', ...
          'Hide blinks in TD plots', 'Truncate TD plot'};
      
fields = {'contour_length', 'pophist_offset', ...
          'cplot_scale_factor', 'contour_bin_size', ...
          'cplot_remove_bleached', 'tdp_max', ...
          'hideBlinksInTDPlots', 'truncate_tdplot' };
currentopt = cellfun( @(x)num2str(opt.(x)), fields, 'UniformOutput',false );

answer = inputdlg(prompt, 'Change display settings', 1, currentopt);
if isempty(answer), return; end  %user hit cancel

% 2. Save new parameter values from user.
for k=1:numel(answer),
    original = opt.(fields{k});
    opt.(fields{k}) = str2double(answer{k});
    
    if islogical(original)
        opt.(fields{k}) = logical( opt.(fields{k}) );
    end
end

opt.fret_axis = -0.1:opt.contour_bin_size:1.2;
opt.contour_bounds = [1 opt.contour_length opt.fretRange];
handles.options = opt;
defaults = opt;  %save as starting values for future calls to makeplots.

% 3. Redraw plots.
plotData(hObject,handles);

end %FUNCTION changeDisplaySettings



function resetSettings2(hObject,~)
% Reset all display settings to their defaults in cascadeConstants.

handles = guidata(hObject);

constants = cascadeConstants;
opt = constants.defaultMakeplotsOptions;
opt.contour_bounds = [1 opt.contour_length opt.fretRange];
handles.options = opt;
defaults = opt;  %save as starting values for future calls to makeplots.

plotData(hObject,handles);

end %FUNCTION resetSettings




end %FUNCTION MAKEPLOTS




function saveFiles2(hObject,~)
% Save plot data to text files for importing and plotting in Origin.
    
handles = guidata(hObject);
base = handles.baseFilenames;

for k=1:numel(base),
    % Save contour FRET histogram
    dlmwrite( [base{k} '_hist.txt'] , handles.cplotdataAll{k}, ' ' );

    % Save display-format (time-binned) contour FRET histogram.
    dlmwrite( [base{k} '_displayhist.txt'], handles.cpdataAll{k}, ' ');

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

histmax = zeros(nFiles,1);
histax = zeros(nFiles,1);  %1D histogram axes
tdax = zeros(nFiles,1);  %tdplot axes

cplotdataAll = cell(nFiles,1); %contour plots, all data
cpdataAll    = cell(nFiles,1); %contour plots, as displayed
shistAll     = cell(nFiles,1); %state occupancy histograms
tdpAll       = cell(nFiles,1); %td plots

% load idealization data from previous calls if possible.
% Warning: this could be too large to fit in memory?
if isfield(handles,'idl'),
    idl = handles.idl;
else
    idl = cell(nFiles,1);
end



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
    
    
    %% ============== DRAW POPULATION CONTOUR HISTOGRAMS ============== 
    
    %---- LOAD OR GENERATE FRET CONTOUR PLOT DATA
    cplotdata = makecplot( fret, options );

    if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
        cplotdata(2:end,2:end) = cplotdata(2:end,2:end).*(N(k)/max(N));
    end
    cplotdataAll{k} = cplotdata;
    
    
    %---- DRAW FRET CONTOUR PLOT ----
    if ~isfield(options,'targetAxes')
        ax = subplot( nrows, nFiles, k );
    else
        ax = options.targetAxes{k,1};
    end
    
    % Draw the contour plot (which may be time-binned)
    % FIXME: we are passing the default options
    cpdataAll{k} = cplot( ax, cplotdata, options.contour_bounds, options );
    
    % Formatting
    title( titles{k}, 'Parent',ax );
    
    if ~options.hideText,
        ap = get(ax,'Position');
        annotation( handles.hFig, 'textbox', [ap(1)+0.9*ap(3) ap(2)+0.93*ap(4) 0.1*ap(3) 0.1*ap(4)], ...
                    'String',sprintf('N=%d', N(k)), 'HorizontalAlignment','right', ...
                    'LineStyle','none' );
    end
    
    if k==1,
         ylabel(ax,'FRET');
         xlabel(ax,'Time (frames)');
    else
        set(ax,'YTickLabel',[])
        ylabel(ax,'');
        xlabel(ax,'');
    end
    ylim( ax, options.fretRange );
    set(ax,'YGrid','on');
    box(ax,'on');
    drawnow;
    
    
    
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
        ax = subplot( nrows, nFiles, nFiles+k );
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
        pophist = sum(histdata,2)/options.contour_length;   %normalization
        
        histmax(k) = max(max( pophist(fretaxis>0.05) ));
        
        bar( ax, fretaxis, pophist );
        

    %---- GENERATE STATE OCCUPANCY HISTOGRAMS
    % ...if dwell-time info is available.
    else
        bins = shist(:,1);
        histdata = shist(:,2:end)*100;
        [~,nStates] = size(histdata);
        
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
    
    
    % Formatting
    hold(ax,'off');
    if k==1,
        ylabel( ax,'Occupancy (%)' );
        xlabel( ax,'FRET' );
    else
       set(ax,'yticklabel',[]);
    end
    set(ax,'YGrid','on');
    xlim( ax,options.fretRange );
    box(ax,'on');
    drawnow;
    
    
    
    %% ========================= DRAW TD PLOTS ========================= 
    
    if ~has_dwt(k) || options.no_tdp
        continue;
    end
    
    
    %---- LOAD TDPLOT DATA ----
    tdp = tdplot(idl{k}, data, options);
    assert( tdp(1,1) ~= 0, 'Invalid TDP' );
    
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
        tdax(k) = subplot( nrows, nFiles, 2*nFiles+k );  
    elseif size(options.targetAxes,2)>2
        tdax(k) = options.targetAxes{k,3};
    else
        continue;
    end
    box(tdax(k),'on');
    tplot( tdax(k), tdp, options );
    
    % Formatting
    if ~options.hideText,
        textOpt = {'HorizontalAlignment','center', 'Parent',tdax(k)};
        text( 0.45,0.9, sprintf('N_t=%.0f',t),            textOpt{:} );
        text( 0.45,0.0, sprintf('t/s=%.2f',t/total_time), textOpt{:} );
    end
    grid(tdax(k),'on');

    if k==1,
        ylabel(tdax(k),'Final FRET');
        xlabel(tdax(k),'Inital FRET');
    else
        set(tdax(k),'yticklabel',[]);
    end
    ylim(tdax(k), options.fretRange );
    xlim(tdax(k), options.fretRange );
    
    
    drawnow;
    
end  %for each file.



%% =================== FINISH UP =================== 

% Save results back to handles object
handles.cplotdataAll = cplotdataAll;
handles.cpdataAll = cpdataAll;
handles.shistAll = shistAll;
handles.tdpAll = tdpAll;
handles.idl = idl;

guidata(hObject,handles);


% Scale histograms to match
if any(histax)
    linkaxes( histax, 'xy' );
    ylim( histax(1), [0 max(histmax)+0.5] );
end

if any(tdax),
    linkaxes( tdax(tdax~=0), 'xy' );
end


set(handles.hFig,'pointer','arrow'); drawnow;

end  % function plotData








