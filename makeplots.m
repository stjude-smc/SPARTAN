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
    
    if isfield(options,'constants')
        constants = options.constants;
    end
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
        titles{i} = strip_path( titles{i} );
    end
    
end


disp('Computing plots, please wait...');


% Determine which files have dwell-time data so that TD plots should be
% made. This determines the number of subplot rows ahead of time.
has_dwt = zeros(nSamples,1);

for i=1:nSamples,
    if exist([baseFilenames{i} '.qub.dwt'],'file'),
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
end


%% ===================== LOOP OVER EACH DATA FILE ====================== 

% If we choose an option that needs to the number of traces in each file,
% we have to load that first before the main loop...
if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
    for i=1:numel(baseFilenames),
        data = loadTraces( dataFilenames{i} );
        N(i) = size( data.fret, 1 );
    end
end

for i=1:numel(baseFilenames),  %for each sample
        
    data_fname  = dataFilenames{i};
    dwt_fname   = [baseFilenames{i} '.qub.dwt'];
    
    
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
    N(i) = size(fret,1);
    
    
    
    %% ============== DRAW POPULATION CONTOUR HISTOGRAMS ============== 
    
    %---- LOAD OR GENERATE FRET CONTOUR PLOT DATA
    if exist('cplotDataArray','var')
        cplotdata = cplotDataArray{i};
        % FIXME: we don't know N here!!
    
    % Generate the contour plot if not available
    else        
        cplotdata = makecplot( fret, options );
        
        if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
            cplotdata(2:end,2:end) = cplotdata(2:end,2:end).*(N(i)/max(N));
        end
        cplotdataAll{i} = cplotdata;
    end
    
    
    %---- DRAW FRET CONTOUR PLOT ----
    if ~isfield(options,'targetAxes')
        ax = subplot( nrows, nSamples, i );
    else
        ax = options.targetAxes{i,1};
        %axes(ax);
    end
    
    % Draw the contour plot (which may be time-binned)
    cpdata = cplot( ax, cplotdata, options.contour_bounds, constants );
    cpdataAll{i} = cpdata;
    
    % Formatting
    title( titles{i}, 'FontSize',16, 'FontWeight','bold', 'Parent',ax );
    
    if ~options.hideText,
        text( 0.90*options.contour_length, 0.94*options.contour_bounds(4), ...
              sprintf('N=%d', N(i)), ...
              'FontWeight','bold', 'FontSize',14, ...
              'HorizontalAlignment','right', 'Parent',ax );
    end
    
    if i==1,
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
    if ~options.no_statehist && has_dwt(i),
        shist = statehist( dwt_fname, fret, options );        
        shistAll{i} = shist;
    end
    
    
    if ~isfield(options,'targetAxes')
        ax = subplot( nrows, nSamples, nSamples+i );
    elseif size(options.targetAxes,2)>1
        ax = options.targetAxes{i,2};
    else
        continue;
    end
    histx(i) = ax;
    cla(ax);  hold(ax,'on');
    set(ax,'ColorOrder',options.colors);
    
    
    %---- GENERATE STATE OCCUPANCY HISTOGRAMS
    % ...if no dwell-time info is available.
    if ~has_dwt(i) || options.no_statehist || isempty(shist)
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
        [nBins,nStates] = size(histdata);
        
        % If requested, remove 0-FRET peak and renormalize
        if options.ignoreState0
            nStates = nStates-1;
            
            histdata = histdata(:,2:end);
            histdata = 100*histdata ./ sum(histdata(:));
        end
        
        % If the option is set, rescale so that plots with only a few
        % molecules show low occupancy in the statehist.
        if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
            histdata = histdata.*(N(i)/max(N));
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
    if i==1,
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
    if ~has_dwt(i) || options.no_tdp
        disp('TDP data missing, skipping');
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
        tdp(2:end,2:end) = tdp(2:end,2:end).*(N(i)/max(N));
    end
    
    tdpAll{i} = tdp;
    
    
    %---- DISPLAY TD PLOT ----
    if ~isfield(options,'targetAxes')
        tdx(i) = subplot( nrows, nSamples, 2*nSamples+i );  
    elseif size(options.targetAxes,1)>2
        tdx(i) = axes(options.targetAxes{i,3});
    else
        disp('no 3rd axis for TD plots...');
        continue;
    end
    box(tdx(i),'on');
    tplot( tdp );
    
    % Formatting
    if ~options.hideText,
        text( 0.43,0.8, sprintf('N_t=%.0f',t), 'FontSize',14, ...
              'FontWeight','bold', 'HorizontalAlignment','center', 'Parent',tdx(i) );
        text( 0.43,0.0, sprintf('t/s=%.1f',t/total_time), 'FontSize',14, ...
              'FontWeight','bold', 'HorizontalAlignment','center', 'Parent',tdx(i) );
    end
    grid(tdx(i),'on');

    if i==1,
        ylabel(tdx(i),'Final FRET');
        xlabel(tdx(i),'Inital FRET');
    else %i>1,
        set(tdx(i),'yticklabel',[]);
    end
    ylim(tdx(i), options.fretRange );
    xlim(tdx(i), options.fretRange );
    
    
    drawnow;
    
end  % for each data file




% Give a warning if using some funky normalization.
if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max
    disp('NOTE: these plots are normalized to the plot with the largest number of traces!!');
end




%% =================== FINISH UP =================== 

% Scale histograms to match
if any(histx)
    linkaxes( histx, 'xy' );
    ylim( histx(1), [0 histmax] );
end

if nrows==3 && numel(tdx>1),
    linkaxes( tdx(tdx~=0), 'xy' );
end



%% =============== ADD GUI CONTROLS ================ 
% These are buttons and things to save the plots or manipulate them without
% calling makeplots again.


% Button to save data
uicontrol( 'Style','pushbutton', 'String','Save files', ...
           'Position',[20 20 80 30], 'Callback',@saveFiles, ...
           'Parent',h1 );


                 
                 
%% ================ GUI CALLBACKS ================ 
% These are defined within the main function scope so we can steal the
% variables (data, filenames, etc).

function saveFiles(hObject,e)
    
    for i=1:numel(baseFilenames),
        base = baseFilenames{i};
    
        % Save contour FRET histogram
        dlmwrite( [base '_hist.txt'] , cplotdataAll{i}, ' ' );
        
        % Save display-format (time-binned) contour FRET histogram.
        dlmwrite( [base '_displayhist.txt'], cpdataAll{i}, ' ');
        
        % State histogram
        if ~isempty(shistAll{i}),
            dlmwrite( [base '.qub_shist.txt'], shistAll{i}, ' ' );
        end
        
        % TD Plot
        if ~isempty(tdpAll{i}),
            dlmwrite( [base '.qub_tdp.txt'], tdpAll{i}, ' ' );
        end
    end
    
end

end %function makeplots
    





%% =============== OTHER FUNCTIONS =============== 


function output = strip_path( filename )

pos = find( filename==filesep );
output = filename(pos(end)+1:end);

end



