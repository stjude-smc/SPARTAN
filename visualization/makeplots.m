function [h1,baseFilenames] = makeplots(varargin)
% function h1 = makeplots(baseFilenames, titles)
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
%

% Created by: Daniel Terry (Scott Blanchard Lab)
% Cascade smFRET Analysis Pipeline 1.3, Copyright (C) 2008 Daniel Terry
% Date Created: Oct 11, 2007

% TODO: add varargin options for stuff




%------- MAKEPLOTS OPTION/PARAMETER VALUES -------

constants = cascadeConstants();
options = constants.defaultMakeplotsOptions;

% Make any changes to the defaults here...

%---------------------------------




%% INITIALIZE & PROCESS FUNCTION PARAMETERS

% Option 1: If no files specified, prompt user for them.
if nargin==0,  %~exist('baseFilenames','var'),
    disp('Select traces files, hit cancel when finished');
    baseFilenames = getFiles([],'Choose a traces file:');

% Option 2: plot data was passed directly (see autotrace.m)
elseif isnumeric(varargin{1})
    cplotDataArray = varargin{1};
    if ~iscell(cplotDataArray)
        cplotDataArray = {cplotDataArray};
    end
    nSamples = length(cplotDataArray);
    baseFilenames  = cell(nSamples,1);
    [baseFilenames{:}] = deal('----');
    options.saveFiles = false;

% Option 3: filenames of data were passed.
else
    baseFilenames = varargin{1};
end

if ~iscell(baseFilenames), baseFilenames={baseFilenames}; end


%% Get optional parameter values

% Titles for each plot set
if nargin>1,
    titles = varargin{2};
    if ~iscell(titles)
        titles = {titles};
    end
end

% Other options (passed as a structure)
if nargin>2,
    options = catstruct( options, varargin{3}, 'sorted' );
end


% Declare persistent variables. If any parameters values are changed in the GUI,
% the new value will override the default until MATLAB is restarted.
persistent persistent_options;

if ~isempty(persistent_options),
    options = catstruct( options, persistent_options, 'sorted' );
end


options.contour_bounds = [1 options.contour_length options.fretRange];


%% Process data to produce plots

nSamples = numel(baseFilenames);

if nSamples == 0,
    disp('No files specified, exiting.');
    return;
end


% Strip extensions from files to get the "base name" from which all other
% filenames can be generated (_tdp.txt _hist.txt, etc)
dataFilenames = baseFilenames;

baseFilenames = strrep( baseFilenames, '.txt', '' );
baseFilenames = strrep( baseFilenames, '.traces', '' );

% Generate plot titles 
if ~exist('titles','var'),
    
    % Remove underscores (subscript)
    titles = strrep(baseFilenames,'_',' ');
    
    % Strip off path, leaving just filename
    for i=1:nSamples,
        [~,titles{i}] = fileparts( titles{i} );
    end
    
end


disp('Computing plots, please wait...');


% Determine which files have dwell-time data so that TD plots should be
% made. This determines the number of subplot rows ahead of time.
has_dwt = zeros(nSamples,1);

for i=1:nSamples,
    if exist([baseFilenames{i} '.qub.dwt'],'file') || ...
       exist([baseFilenames{i} '.dwt'],'file'),
        has_dwt(i) = true;
    end
end %for each sample

if any(has_dwt) && ~options.no_tdp, 
    nrows = 3;
else
    nrows = 2;
end


histmax = 0.05;
histx = zeros(nSamples,1);  %1D histogram axes
tdx = zeros(nSamples,1);  %tdplot axes
N = zeros(nSamples,1); %number of molecules, each sample

cplotdataAll = cell(nSamples,1); %contour plots, all data
cpdataAll    = cell(nSamples,1); %contour plots, as displayed
shistAll     = cell(nSamples,1); %state occupancy histograms
tdpAll       = cell(nSamples,1); %td plots

if ~isfield(options,'targetAxes')
    h1 = figure();
else                                                    %MJ
    h1 = get(options.targetAxes{1,1},'parent');         %MJ
end


%% ===================== LOOP OVER EACH DATA FILE ====================== 

% If we choose an option that needs the number of traces in each file,
% we have to load that first before the main loop...
if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
    for i=1:numel(baseFilenames),
        N(i) = sizeTraces( dataFilenames{i} );
    end
end

% See the end of the file. This loads the data, calculates plots, and draws
% them.
plotData();


% Give a warning if using some funky normalization.
if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max
    disp('NOTE: these plots are normalized to the plot with the largest number of traces!!');
end




%% =============== ADD GUI CONTROLS ================ 
% These are buttons and things to save the plots or manipulate them without
% calling makeplots again.
% Position from LL corner is defined as [left bottom width height].
% FIXME: change to normalized units so the buttons are scaled?

uicontrol( 'Style','pushbutton', 'String','Save files', ...
           'Position',[20 20 80 30], 'Callback',@saveFiles, ...
           'Parent',h1 );

uicontrol( 'Style','pushbutton', 'String','Change settings', ...
           'Position',[110 20 130 30], 'Callback',@changeDisplaySettings, ...
           'Parent',h1 );
       
uicontrol( 'Style','pushbutton', 'String','Reset settings', ...
           'Position',[250 20 120 30], 'Callback',@resetSettings, ...
           'Parent',h1 );

       
       
       
                 
                 
%% ================ GUI CALLBACKS ================ 
% These are defined within the main function scope so we can steal the
% variables (data, filenames, etc).
% FIXME: consider having a single dialog for all major settings.

function saveFiles(~,~)
% Save plot data to text files for importing and plotting in Origin.
    
    for k=1:numel(baseFilenames),
        base = baseFilenames{k};
    
        % Save contour FRET histogram
        dlmwrite( [base '_hist.txt'] , cplotdataAll{k}, ' ' );
        
        % Save display-format (time-binned) contour FRET histogram.
        dlmwrite( [base '_displayhist.txt'], cpdataAll{k}, ' ');
        
        % State histogram
        if ~isempty(shistAll{k}),
            dlmwrite( [base '.qub_shist.txt'], shistAll{k}, ' ' );
        end
        
        % TD Plot
        if ~isempty(tdpAll{k}),
            dlmwrite( [base '.qub_tdp.txt'], tdpAll{k}, ' ' );
        end
    end
    
end %FUNCTION saveFiles



function changeDisplaySettings(~,~)
% Changes how much (how many frames) of the movie to show in plots.
% This is equivalent to changing pophist_sumlen in cascadeConstants.


% 1. Get the new value from the user.
answer = inputdlg( {'Contour length (frames):','Contour offset (frames):', ...
                    'Contour scaling factor:','TD plot scaling factor:' }, ...
         'Change display settings', 1, ...
        { num2str(options.contour_length),     num2str(options.pophist_offset), ...
          num2str(options.cplot_scale_factor), num2str(options.tdp_max) }  );

if isempty(answer), return; end  %user hit cancel


% 2. Save new parameter values from user.
persistent_options.contour_length = str2double( answer{1} );
persistent_options.pophist_offset = str2double( answer{2} );
persistent_options.cplot_scale_factor = str2double( answer{3} );
persistent_options.tdp_max            = str2double( answer{4} );

options = catstruct( options, persistent_options, 'sorted' );
options.contour_bounds = [1 options.contour_length options.fretRange];


% 3. Redraw plots.
plotData();


end %FUNCTION changeDisplaySettings



function resetSettings(~,~)
% Reset all display settings to their defaults in cascadeConstants. Any
% persistent settings will be overwritten

persistent_options = [];

constants = cascadeConstants;
options = constants.defaultMakeplotsOptions;
options.contour_bounds = [1 options.contour_length options.fretRange];

plotData();

end %FUNCTION changeDisplaySettings






%% =======================================================================

%% ===================== LOOP OVER EACH DATA FILE ====================== 

function plotData()
% This is the function that actually loads the data, calculates the plots,
% and displays them.

for k=1:numel(baseFilenames),  %for each sample
        
    data_fname  = dataFilenames{k};
    dwt_fname   = [baseFilenames{k} '.qub.dwt'];
    
    if has_dwt(k) && ~exist(dwt_fname,'file'),
        dwt_fname = [baseFilenames{k} '.dwt'];
    end
    
    % Load FRET data
    if ~exist(data_fname,'file')
        disp('FRET Data missing, skipping.'); continue;
    end
    
    data = loadTraces( data_fname );
    
    % The data may have multiple FRET signals to analyze. Ask which to use.
    if data.isChannel('fret2'),
        a = questdlg('This data has multiple FRET channels. Which should be used?', ...
                     'Select FRET channel to use','fret','fret2','Cancel', 'fret');
        if strcmp(a,'Cancel'), return; end
        fret = data.(a);
    else
        fret = data.fret;
    end
    
    clear data;
    N(k) = size(fret,1);
    
    
    
    %% ============== DRAW POPULATION CONTOUR HISTOGRAMS ============== 
    
    %---- LOAD OR GENERATE FRET CONTOUR PLOT DATA
    if exist('cplotDataArray','var')
        cplotdata = cplotDataArray{k};
        % FIXME: we don't know N here!!
    
    % Generate the contour plot if not available
    else
        cplotdata = makecplot( fret, options );
        
        if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
            cplotdata(2:end,2:end) = cplotdata(2:end,2:end).*(N(k)/max(N));
        end
        cplotdataAll{k} = cplotdata;
    end
    
    
    %---- DRAW FRET CONTOUR PLOT ----
    if ~isfield(options,'targetAxes')
        ax = subplot( nrows, nSamples, k );
    else
        ax = options.targetAxes{k,1};
        %axes(ax);
    end
    
    % Draw the contour plot (which may be time-binned)
    % FIXME: we are passing the default options
    cpdata = cplot( ax, cplotdata, options.contour_bounds, options );
    cpdataAll{k} = cpdata;
    
    % Formatting
    title( titles{k}, 'FontSize',16, 'FontWeight','bold', 'Parent',ax );
    
    if ~options.hideText,
        text( 0.90*options.contour_length, 0.94*options.contour_bounds(4), ...
              sprintf('N=%d', N(k)), ...
              'FontWeight','bold', 'FontSize',14, ...
              'HorizontalAlignment','right', 'Parent',ax );
    end
    
    if k==1,
         ylabel(ax,'FRET');
         xlabel(ax,'Time (frames)');
    else %i>1
        set(ax,'YTickLabel',[])
        ylabel(ax,'');
        xlabel(ax,'');
    end
    ylim( ax, options.fretRange );
    set(ax,'YGrid','on');
    box(ax,'on');
    
    
    
      
    %% ================ DRAW STATE OCCUPANCY HISTOGRAMS ================ 
   
    %---- GENERATE STATE OCCUPANCY HISTOGRAMS
    if ~options.no_statehist && has_dwt(k),
        shist = statehist( dwt_fname, fret, options );        
        shistAll{k} = shist;
    end
    
    
    if ~isfield(options,'targetAxes')
        ax = subplot( nrows, nSamples, nSamples+k );
    elseif size(options.targetAxes,2)>1
        ax = options.targetAxes{k,2};
    else
        continue;
    end
    histx(k) = ax;
    cla(ax);  hold(ax,'on');
    set(ax,'ColorOrder',options.colors);
    
    
    %---- GENERATE STATE OCCUPANCY HISTOGRAMS
    % ...if no dwell-time info is available.
    if ~has_dwt(k) || options.no_statehist || isempty(shist)
        fretaxis = cplotdata(2:end,1);      
        histdata = cplotdata(2:end,2:options.contour_length+1)*100;
        pophist = sum(histdata,2)/options.contour_length;   %normalization
        
        max_bar = max(max( pophist(fretaxis>0.05) ));
        
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
        
        max_bar = max( totalHist(bins>0.05) );
    end
    
    
    % Histogram height axis adjustment
    histmax = max( histmax, max_bar+0.5 );
    
    
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
    
    
    
    %% ========================= DRAW TD PLOTS ========================= 
    
    % Make sure TD plot data exists
    if ~has_dwt(k) || options.no_tdp
        %disp('TDP data missing, skipping');
        continue;
    end
    
    
    %---- GENERATE TD PLOT HISTOGRAMS
    tdp = tdplot(dwt_fname,fret,options);
 
    
    %---- LOAD TDPLOT DATA ----
    
    % If total time (normalization factor) is saved in TD plot, use it to
    % get raw number of transitions
    if tdp(1,1) ~= 0,
        total_time = tdp(1,1);
        t = sum(sum( tdp(2:end,2:end)*total_time ));
        
    % Otherwise, calculate it from the idealization file (DWT)
    else
        if exist(dwt_fname,'file')~=2,
            disp('DWT needed but missing, skipping');
            continue;
        end
        
        [transmat,total_time] = CountEvents(dwt);
        t = sum(sum( transmat ));
    end
    
        
    % If the option is set, rescale so that plots with only a few
    % molecules show low occupancy in the statehist.
    if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
        tdp(2:end,2:end) = tdp(2:end,2:end).*(N(k)/max(N));
    end
    
    tdpAll{k} = tdp;
    
    
    %---- DISPLAY TD PLOT ----
    if ~isfield(options,'targetAxes')
        tdx(k) = subplot( nrows, nSamples, 2*nSamples+k );  
    elseif size(options.targetAxes,2)>2                 %MJ
        tdx(k) = options.targetAxes{k,3};
%         options.ax = options.targetAxes{k,3};           %MJ
    else
        disp('no 3rd axis for TD plots...');
        continue;
    end
    box(tdx(k),'on');
    tplot( tdx(k), tdp, options );
    
    % Formatting
    if ~options.hideText,
        text( 0.43,0.8, sprintf('N_t=%.0f',t), 'FontSize',14, ...
              'FontWeight','bold', 'HorizontalAlignment','center', 'Parent',tdx(k) );
        text( 0.43,0.0, sprintf('t/s=%.2f',t/total_time), 'FontSize',14, ...
              'FontWeight','bold', 'HorizontalAlignment','center', 'Parent',tdx(k) );
    end
    grid(tdx(k),'on');

    if k==1,
        ylabel(tdx(k),'Final FRET');
        xlabel(tdx(k),'Inital FRET');
    else %i>1,
        set(tdx(k),'yticklabel',[]);
    end
    ylim(tdx(k), options.fretRange );
    xlim(tdx(k), options.fretRange );
    
    
    drawnow;
    
end  %for each file.


%% =================== FINISH UP =================== 

% Scale histograms to match
if any(histx)
    linkaxes( histx, 'xy' );
    ylim( histx(1), [0 histmax] );
end

if nrows==3 && numel(tdx>1),
    linkaxes( tdx(tdx~=0), 'xy' );
end

    

end  % function plotData











end %function makeplots.

