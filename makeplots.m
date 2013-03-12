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
%   OPTIONS is a structure with any of the following fields:
%     .pophistSumlen = number of frames to plot for histograms
%     .contourBinSize = 
%     .cplotMax       = threshold for highest level in contour plots
%     .noStatehist    = 1: do not make/show state histograms
%     .noTDP          = 1: do not make/show transition-density plots
%     .ignoreState0   = 1: do not show zero-FRET state in state histograms
%     .hideText       = 1: do not show text on histograms
%     .targetAxes     = cell array of handles to axes in which to plot
%   Default values are shown below under "default parameter values".

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


% EXTENSIONS for files used:
dwt_ext  = '.qub.dwt';       % QuB Idealization file
hist_ext = '_hist.txt';      % 1D population histogram
tdp_ext  = '.qub_tdp.txt';   % TD plot
shist_ext= '.qub_shist.txt';     % state occupancy histogram





%% =================== CREATE TD & STATEHIST PLOTS =================== 

disp('Computing plots, please wait...');

any_tdps = false;  % will we be drawing any TD plots?

for i=1:nSamples,
    
    data_fname  = dataFilenames{i};
    dwt_fname   = [baseFilenames{i} dwt_ext];
    
    % Make sure files exist
    any_tdps = 0;
    
    if ~exist(dwt_fname,'file')
        disp('Idealization data missing'); continue;
    elseif ~exist(data_fname,'file')
        disp('FRET Data missing, skipping.'); continue;
    else
        any_tdps = true;
    end
    
    %---- GENERATE TD PLOT HISTOGRAMS
    if ~options.no_tdp,
        tdplot(dwt_fname,data_fname,options);  % save file: '_tdp.txt' extension
    end
    
    
    %---- GENERATE STATE OCCUPANCY HISTOGRAMS
    if ~options.no_statehist,
        statehist( dwt_fname, data_fname, options );
    end
end

if options.no_tdp, any_tdps = 0; end

if any_tdps, 
    nrows = 3;  %only draw population histograms
else
    nrows = 2;  %also draw TD plots
end


%% =================== DRAW POPULATION HISTOGRAMS =================== 

histmax = 0.05;
histx = zeros(nSamples,1);


if ~isfield(options,'targetAxes')
    h1 = figure();
end
    
% set(gcf,'DefaultAxesColorOrder',colors);

N = zeros(nSamples,1); %number of molecules, each sample

for i=1:numel(baseFilenames),  %for each sample
    
    hist_fname = [baseFilenames{i} hist_ext];
    displayhist_fname = [baseFilenames{i} '_displayhist.txt'];
    shist_fname   = [baseFilenames{i} shist_ext];
    data_fname = dataFilenames{i};
    
    %---- LOAD OR GENERATE FRET CONTOUR PLOT DATA
    
    if exist('cplotDataArray','var')
        cplotdata = cplotDataArray{i};
        % FIXME: we don't know N here!!
    
    % Generate the contour plot if not available
    else
        % Make sure FRET data exists
        if ~exist(data_fname,'file')
            disp('Traces file missing, skipping');
            continue;
        end
        
        disp('Generating colorplot...');
        
        data = loadTraces(data_fname);
        N(i) = size(data.donor,1);
        
        cplotdata = makecplot( data.fret, options );
        dlmwrite(hist_fname,cplotdata,' ');
    end
    
    
    %---- DRAW FRET CONTOUR PLOT ----
    
    if ~isfield(options,'targetAxes')
        ax = subplot( nrows, nSamples, i );
    else
        ax = options.targetAxes{i,1};
        axes(ax);
    end
    
    % Draw the contour plot (which may be time-binned)
    cpdata = cplot( ax, cplotdata, options.contour_bounds, constants );
    
    % Save display-format (time-binned) histogram.
    if exist(data_fname,'file'),
        dlmwrite(displayhist_fname,cpdata,' ');
    end
    
    % Formatting
    title( titles{i}, 'FontSize',16, 'FontWeight','bold', 'Parent',ax );
    
    if ~options.hideText && N(i)>0,
        axes(ax);
        text( 0.93*options.contour_length, 0.9*options.contour_bounds(4), ...
              sprintf('N=%d', N(i)), ...
              'FontWeight','bold', 'FontSize',14, ...
              'HorizontalAlignment','right', 'Parent',ax );
    end
    
    if i==1,
         ylabel(ax,'FRET');
         xlabel(ax,'Time (frames)');
%          set(gca,'ytick',0:0.2:1.0);
    else %i>1
        set(ax,'YTickLabel',[])
        ylabel(ax,'');
        xlabel(ax,'');
    end
    ylim( ax, options.fretRange );
    set(ax,'YGrid','on');
    box(ax,'on');
    
      
    
    %---- DRAW STATE OCCUPANCY HISTOGRAM ----
   
    if ~isfield(options,'targetAxes')
        ax = subplot( nrows, nSamples, nSamples+i );
    elseif size(options.targetAxes,2)>1
        ax = options.targetAxes{i,2};
        axes(ax);
    else
        continue;
    end
    histx(i) = ax;
    cla(ax);  hold(ax,'on');
    set(ax,'ColorOrder',options.colors);
    
    % If idealization data missing, plot 1D population histogram
    if ~exist(shist_fname,'file') || options.no_statehist
        fretaxis = cplotdata(2:end,1);      
        histdata = cplotdata(2:end,2:options.contour_length+1)*100;
        pophist = sum(histdata,2)/options.contour_length;   %normalization
        
        max_bar = max(max( pophist(fretaxis>0.05) ));
        
        bar( ax, fretaxis, pophist );


    % Otherwise, plot state occupancy histogram
    else
        shist = load( shist_fname );
        bins = shist(:,1);
        histdata = shist(:,2:end)*100;
        [nBins,nStates] = size(histdata);
        
        % If requested, remove 0-FRET peak and renormalize
        if options.ignoreState0
            nStates = nStates-1;
            
            histdata = histdata(:,2:end);
            histdata = 100*histdata ./ sum(histdata(:));
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
    
end


% Scale histograms to match
if exist('histx','var')
    linkaxes( histx, 'xy' );
    ylim( histx(1), [0 histmax] );
end

drawnow;




%% =================== DRAW TD PLOTS =================== 

if nrows < 3,  return;  end
tdx = zeros(nSamples,1);

for i=1:nSamples,  %for each sample
    
    tdp_fname = [baseFilenames{i} tdp_ext];
    dwt_fname = [baseFilenames{i} dwt_ext];
    
    
    % Make sure TD plot data exists
    if ~exist(tdp_fname,'file')
        disp('TDP data missing, skipping');
        continue;
    end
    
    
    %---- LOAD TDPLOT DATA ----
    tdp = load(tdp_fname);
    
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
%         set(gca,'ytick',1:0.2:0.9);
    else %i>1,
        set(tdx(i),'yticklabel',[]);
    end
    ylim(tdx(i), options.fretRange );
    xlim(tdx(i), options.fretRange );
end

if numel(tdx>1) && all(tdx~=0),
    linkaxes( tdx, 'xy' );
end
drawnow;

end %function makeplots
    



%% =============== OTHER FUNCTIONS =============== 


function output = strip_path( filename )

pos = find( filename==filesep );
output = filename(pos(end)+1:end);

end



