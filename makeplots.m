function [h1,samples] = makeplots(varargin)
% function h1 = makeplots(samples, titles)
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




%------- DEFAULT PARAMETERS VALUES -------

constants = cascadeConstants();

% Range of FRET values to display in countor and TD plots
options.fretRange = [-0.1 0.99];

% Length (in frames) of population histograms
options.pophist_sumlen = constants.contour_length;

% axis bounds for contour plot and 1D population histogram
options.contour_bin_size = 0.03;

% colors for statehist, in order
options.colors = [ 0 0 0 ; 0 0.5 0   ; 1 0 0 ; ...
               0    0.7500    0.7500 ; ...
          0.7500         0    0.7500 ; ...
          0.7500    0.7500         0 ; ...
          0.2500    0.2500    0.2500 ];

% OPTIONS
options.force_remake_tdplot = false;  %regenerate statehist and tdplots
options.no_statehist = true;  % do not use state occupancy histograms
options.no_tdp       = true;  % do not use TD plots
options.ignoreState0 = true;   % do not include the first (lowest FRET) state in
                       % state occupancy histograms
options.hideText     = false;  % don't display N=, t/s, etc on plots
                       
if options.ignoreState0,
   options.colors = options.colors(2:end,:); 
end

%---------------------------------

















%% INITIALIZE & PROCESS FUNCTION PARAMETERS

% Option 1: If no files specified, prompt user for them.
if nargin==0,  %~exist('samples','var'),
    disp('Select traces files, hit cancel when finished');
    samples = getFiles('*.txt','Choose a traces file:');

% Option 2: plot data was passed directly (see autotrace.m)
elseif isnumeric(varargin{1})
    cplotDataArray = varargin{1};
    if ~iscell(cplotDataArray)
        cplotDataArray = {cplotDataArray};
    end
    nSamples = length(cplotDataArray);
    samples  = cell(nSamples,1);
    [samples{:}] = deal('----');

% Option 3: filenames of data were passed.
else
    samples = varargin{1};
end


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

options.contour_bounds = [1 options.pophist_sumlen options.fretRange];


%% Process data to produce plots

nSamples = numel(samples);

if nSamples == 0,
    disp('No files specified, exiting.');
    return;
end


% Strip extensions from files to get the "base name" from which all other
% filenames can be generated (_tdp.txt _hist.txt, etc)
samples = strrep( samples, '.txt', '' );
samples = strrep( samples, '_tdp', '' );
samples = strrep( samples, '.qub', '' );
samples = strrep( samples, '_hist', '' );


% Generate plot titles 
if ~exist('titles','var'),
    
    % Remove underscores (subscript)
    titles = strrep(samples,'_',' ');
    
    % Strip off path, leaving just filename
    for i=1:nSamples,
        titles{i} = strip_path( titles{i} );
    end
    
end


% EXTENSIONS for files used:
dwt_ext  = '.qub.dwt';       % QuB Idealization file
% data_ext = '.qub.txt';       % forQuB raw data (not used)
data_ext  = '.txt';           % raw traces data
hist_ext = '_hist.txt';      % 1D population histogram
tdp_ext  = '.qub_tdp.txt';   % TD plot
shist_ext= '_shist.txt';     % state occupancy histogram





%% =================== CREATE TD & STATEHIST PLOTS =================== 

disp('Computing plots, please wait...');

any_tdps = false;  % will we be drawing any TD plots?

for i=1:nSamples,
    
    data_fname  = [samples{i} data_ext];
    dwt_fname   = [samples{i} dwt_ext];
    tdp_fname   = [samples{i} tdp_ext];
    shist_fname = [samples{i} shist_ext];
    
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
    
    if options.no_tdp, continue; end
    
    % Speed: create TD plot hist only if new or if data files have changed
    datadir = dir(data_fname);
    dwtdir  = dir(dwt_fname);
    tdpdir =  dir(tdp_fname);
    
    data_date = max( datadir.datenum, dwtdir.datenum );
    
    if numel(tdpdir)==0 || data_date>tdpdir.datenum || options.force_remake_tdplot
        disp('New data detected: generating TD plot hist...');
        tdplot(dwt_fname,data_fname);  % save file: '_tdp.txt' extension
    end
    
    
    %---- GENERATE STATE OCCUPANCY HISTOGRAMS
    if options.no_statehist, continue; end
    
    shdir = dir(shist_fname);
    
    if numel(shdir)==0 || data_date>shdir.datenum || options.force_remake_tdplot
        fretaxis = contour_bounds(3):contour_bin_size:contour_bounds(4);
        
        disp('New data detected: generating state hist...');
        shist = statehist( dwt_fname, data_fname, fretaxis );
        save( shist_fname, 'shist', '-ASCII' );
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


if ~isfield(options,'targetAxes')
    h1 = figure();
end
    
% set(gcf,'DefaultAxesColorOrder',colors);

N = zeros(1,nSamples); %number of molecules, each sample

for i=1:numel(samples),  %for each sample
    
    hist_filename = [samples{i} hist_ext];
    shist_fname   = [samples{i} shist_ext];
    data_fname    = [samples{i} data_ext];
    
    %---- LOAD OR GENERATE FRET CONTOUR PLOT DATA
    
    if exist('cplotDataArray','var')
        cplotdata = cplotDataArray{i};
    
    % Generate the contour plot if not available
    elseif ~exist(hist_filename,'file') || fileIsNewer(data_fname,hist_filename)
        % Make sure FRET data exists
        filename = [samples{i} data_ext];
        if ~exist(filename,'file')
            disp('Traces file missing, skipping');
            continue;
        end
        
        disp('Generating colorplot...');
        cplotdata = makecplot( filename, options.contour_bin_size );
    else
        cplotdata = load(hist_filename);
    end
    
    N(i) = sum( cplotdata(2:end,2) );  %number of traces
    
    
    
    %---- DRAW FRET CONTOUR PLOT ----
    
    if ~isfield(options,'targetAxes')
        ax = subplot( nrows, nSamples, i );
    else
        ax = options.targetAxes{i,1};
        axes(ax);
    end
    
    % Draw the contour plot
    cplot( ax, cplotdata, options.contour_bounds, constants );
    
    % Formatting
    title( titles{i}, 'FontSize',16, 'FontWeight','bold', 'Parent',ax );
    
    if ~options.hideText,
        axes(ax);
        text( 0.93*options.pophist_sumlen, 0.9*options.contour_bounds(4), ...
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
        continue
    end
    histx(i) = ax;
    cla(ax);  hold(ax,'on');
    set(ax,'ColorOrder',options.colors);
    
    % If idealization data missing, plot 1D population histogram
    if ~exist(shist_fname,'file') || options.no_statehist
        fretaxis  = cplotdata(2:end,1);      
        histdata = cplotdata(2:end,2:options.pophist_sumlen+1)*100;
        pophist = sum(histdata,2)/N(i)/options.pophist_sumlen;   %normalization
        
        max_bar = max(max( pophist(fretaxis>0.05) ));
        
        bar( ax, fretaxis, pophist );


    % Otherwise, plot state occupancy histogram
    % FIXME: remove dark state (0) from curve -- at least as an option
    % -- and renormalize without it.
    else
        shist = load( shist_fname );
        bins = shist(:,1);
        histdata = shist(:,2:end)*100;
        [nBins,nStates] = size(histdata);
        
        % If requested, remove 0-FRET peak and renormalize
        if ignoreState0
            nStates = nStates-1;
            
            histdata = histdata(:,2:end);
            histdata = 100*histdata ./ sum(histdata(:));
        end
        
        minBin = find( bins>0.125, 1, 'first' );
        combined = max(histdata');
        max_bar = max(combined(minBin:end));

        % Draw translucent, filled area underneath curves
        for j=1:nStates
            patch( bins, histdata(:,j), options.colors(j,:), ...
                    'EdgeColor','none','FaceAlpha',0.25, 'Parent',ax );
        end
        
        % Draw occupancy histograms as solid lines
        plot( ax, bins, histdata, 'LineWidth', 1.5 );
        
        % Calculate percent time spent in each non-zero FRET state
%         occupancy = sum( histdata(:,2:end) ); %ignore zero-FRET state
%         occupancy = occupancy/sum(occupancy); %normalize
        
        % Draw the percent times
        
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
    ylim( ax,options.fretRange );
    box(ax,'on');
    
end


% Scale histograms to match
if exist('histx','var')
    axis( histx, [options.contour_bounds(3:4) 0 histmax] );
    linkaxes( histx, 'xy' );
end

drawnow;




%% =================== DRAW TD PLOTS =================== 

if nrows < 3,  return;  end
tdx = [];

for i=1:nSamples,  %for each sample
    
    tdp_fname = [samples{i} tdp_ext];
    dwt_fname = [samples{i} dwt_ext];
    
    
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
    if ~hideText,
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
end

if numel(tdx>1) && all(tdx~=0),
    linkaxes( tdx, 'xy' );
end
drawnow;

end %function makeplots
    



%% =============== FCN TO MAKE CONTOUR PLOTS =============== 

function frethist = makecplot( data_filename, contour_bin_size)
% MAKECPLOT   creates _hist.txt FRET histogram file

% Load data
[d,a,fret] = loadTraces( data_filename );
[Nmol len] = size(fret);

% Axes for histogram includes all possible data
time_axis = 1:len;
fret_axis = -0.1:contour_bin_size:1.0;

% Initialize histogram array, setting the time step in the first row,
% and the FRET bins in the first column. This is done for import into
% Origin.
frethist = zeros( length(fret_axis)+1, length(time_axis)+1 );
frethist(1,2:end) = time_axis;
frethist(2:end,1) = fret_axis';

frethist(2:end,2:end) = hist( fret, fret_axis  );


% Save plots to file
histfile=strrep(data_filename,'.txt','_hist.txt');
dlmwrite(histfile,frethist,' ');


end %function makecplot



function output = strip_path( filename )

pos = find( filename==filesep );
output = filename(pos(end)+1:end);

end



function answer = fileIsNewer( filename1, filename2 )
% FILEISNEWER  compare two file modification dates
% 
%   BOOL = FILEISNEWER( FILE1, FILE2 )
%   Full path must be provided and files must exist.

dir1 = dir(filename1);
dir2 = dir(filename2);

answer = (dir1.datenum > dir2.datenum);

end






