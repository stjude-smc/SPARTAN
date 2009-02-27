function [models,fits] = frethistComparison2(files,titles)
% FRETHISTCOMPARISON  Overlays multiple 1D FRET histograms for comparison
%
%   Prompts user for locations of traces files to load.  Creates 1D
%   FRET histograms.  Plots each seperately as stair plots, with
%   different colors and a legend.
%   

% OPTIONS:
forceMakePlots = 1;  %ignore existing files


% EXTENSIONS for files used:
hist_ext = '_hist.txt';      % 1D population histogram

% constants:
contour_bin_size = 0.025;
pophist_sumlen = 50;

colors = [ 1     0    0    ; ...
           0     0    1    ; ...
          0   0.6    0     ; ...
          0   0 0];



% Prompt user for filenames if not supplied
if nargin < 1,
    disp('Select traces files, hit cancel when finished');
    files = getFiles('*.txt','Choose a traces file:');
end

nFiles = numel(files);

if nFiles == 0,
    disp('No files specified, exiting.');
    return;
end



[nrow,ncol] = size(files);


% f1 = figure;
f1 = gcf;
set(f1,'defaultaxesfontsize',16);
set(f1,'defaulttextfontsize',18);

% set( gcf, 'Position', [6 563 1590 341] );
% set( gcf, 'Position', [ 1977         365         405         757 ] );

% Plot all of the histograms together for easy comparison
models = struct([]);
fits = cell(nFiles,1);
maxY = 0;

for i=1:nFiles
    
    base = strrep( files{i}, '.txt', '' );
    data_fname    = files{i};
    hist_filename = [base hist_ext];
    
    % Generate the contour plot if not available
    if forceMakePlots || ~exist(hist_filename,'file') || fileIsNewer(data_fname,hist_filename)
        makecplot( data_fname, contour_bin_size );
    end
    
    % Load histogram data
    cplotdata = load(hist_filename);
    N = sum( cplotdata(2:end,2) );  %number of traces
    
    % Plot histogram
    fretaxis = cplotdata(2:end,1);      
    histdata = cplotdata(2:end,2:pophist_sumlen+1)*100;
    pophist  = sum(histdata,2)/N/pophist_sumlen;   %normalization
    maxY = max( maxY, max(pophist(fretaxis>0.1)) );

%     plot( fretaxis+i*0.003, pophist, 'Color',colors(i,:), 'LineWidth',2 );
%     hold on;
end

% hold off;
% xlim( [-0.1 1] );
% set(gca,'xtick',0:0.2:1);
% ylimits = get(gca,'ylim');
ylimits = [0 1.05*maxY];
% ylim( ylimits );
% xlabel('FRET Efficiency');

grid = 1:(nrow*ncol);
grid = reshape(grid,ncol,nrow)';

% Plot FRET histograms individually with gaussian fits
for i=1:nFiles
    col = floor((i-1)/nrow)+1;
    row = mod(i-1,nrow)+1;
    
    base = strrep( files{i}, '.txt', '' );
    hist_filename = [base hist_ext];
    
    % Load histogram data
    cplotdata = load(hist_filename);
    N = sum( cplotdata(2:end,2) );  %number of traces
    
    % Plot histogram
    fretaxis = cplotdata(2:end,1);      
    histdata = cplotdata(2:end,2:pophist_sumlen+1)*100;
    pophist  = sum(histdata,2)/N/pophist_sumlen;   %normalization

    subplot(nrow,ncol,grid(row,col));
    grey = [0.9 0.9 0.9];
    bar( fretaxis, pophist, 1, 'FaceColor',grey, 'EdgeColor','k' );
    hold on;
    %stairs( fretaxis, pophist, 'k' );
    
    % Fit the data to a sum of gaussians
    nComp = 3;
%     if i<=6, nComp = 2;  end
    [result,model] = fitPeaks( fretaxis, pophist, nComp );
    fits{i} = model;
    models = [models result];
    
    % Add the fit to the plot    
    for j=1:nComp,
        a = model.( ['a' num2str(j)] );
        b = model.( ['b' num2str(j)] );
        c = model.( ['c' num2str(j)] );
        h = plot( fretaxis, a.*exp(-((fretaxis-b)./c).^2), ...
                  'k--', 'LineWidth',1.5 );
    end
    
    h = plot(model);
    set(h, 'LineWidth',2 );
    set(h, 'Color',colors(row,:) );
    
    legend off;
    xlim( [-0.1 1] );
    ylim( ylimits );
    set(gca,'xtick',0:0.2:1);
    if row~=nrow,
        set(gca,'xticklabel',{''});
    end
    hold off;
    
    if row==nrow,
        xlabel( 'FRET Efficiency' );
    else
        xlabel('');
    end
    
    if col<2,
        ylabel( 'Occupancy (%)' );
    else
        ylabel('');
        set(gca,'yticklabel',{});
    end
    
    ylimits = get(gca,'ylim');
    ylimits = [0 1.05*maxY];
    ylim( ylimits );
end




if nargin>1,
    legend( titles );
end

end




%% FUNCTIONS

function frethist = makecplot( data_filename, contour_bin_size )
% MAKECPLOT   creates _hist.txt FRET histogram file
% taken from makeplots

% Load data
[d,a,fret] = LoadTraces( data_filename );
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

fret( fret==0 ) = NaN;

frethist(2:end,2:end) = hist( fret, fret_axis  );


% Save plots to file
histfile=strrep(data_filename,'.txt','_hist.txt');
dlmwrite(histfile,frethist,' ');


end %function makecplot



function answer = fileIsNewer( filename1, filename2 )
% FILEISNEWER  compare two file modification dates
% 
%   BOOL = FILEISNEWER( FILE1, FILE2 )
%   Full path must be provided and files must exist.

dir1 = dir(filename1);
dir2 = dir(filename2);

answer = (dir1.datenum > dir2.datenum);

end



function [results,model] = fitPeaks( fret_axis, histdata, nComp )
% 

% if nargin<3,
%     nComp = 3;
% end


% Fit to a sum of gaussians
% l = [  0 -0.01 0   0 0.1 0   0 0.1 0 ];
% l = [  0 -0.01 0   0 0.1 0   0 0.7 0 ];
% l = [  0 0.1 0   0 0.7 0 ];
% u = [  50 0.15 0.15  repmat([25 1 0.12],[1 nComp-1])  ];
u = repmat([25 1 0.13],[1 nComp]);
% s = [  1 0.01 0.05  repmat([1 0.5 0.1],[1 nComp-1])  ];
% s = [  1 0.01 0.1  1 0.55 0.1  1 0.75 0.1  ];
s = [ 1 0.3 0.07  1 0.6 0.07   1 0.8 0.07 ];

% s = s( 1:(nComp*3) );
% u = u( 1:(nComp*3) );

model = fit( fret_axis, histdata, ['gauss' num2str(nComp)], ...
             'Upper', u, 'StartPoint', s );  %'Lower',l

% Save the fit in results, ignoring first component (dark state)
% normFactor = sum( [model.a1, model.a2] );
% results.occupancy = [model.a1 model.a2]./normFactor;
% results.avgFret = [model.b1 model.b2];
% results.stdFret = [model.c1 model.c2];
% 
% results.avgDist = ((1./results.avgFret)-1).^(1/6).*50-16;
%     
% results.area = results.occupancy .* results.stdFret;
results = struct( [] );

end %function


