function f = frethistComparison(titles, files)
% FRETHISTCOMPARISON  Overlays multiple 1D FRET histograms for comparison
%
%   Prompts user for locations of traces files to load.  Creates 1D
%   FRET histograms.  Plots each seperately as stair plots, with
%   different colors and a legend.
%   


% OPTIONS:
forceMakePlots = 1;  %ignore existing files


% EXTENSIONS for files used:
dwt_ext  = '.qub.dwt';       % QuB Idealization file
data_ext = '.qub.txt';       % forQuB raw data (not used)
raw_ext  = '.txt';           % raw traces data
hist_ext = '_hist.txt';      % 1D population histogram
tdp_ext  = '.qub_tdp.txt';   % TD plot
shist_ext= '_shist.txt';     % state occupancy histogram

% constants:
contour_bin_size = 0.025;
pophist_sumlen = 50;

colors = [ 0         0.7500    0.7500 ; ...
           0.7500         0    0.7500 ; ...
           0.7500    0.7500         0 ; ...
           0.2500    0.2500    0.2500 ];



% Prompt user for filenames if not supplied
if nargin < 2,
    disp('Select traces files, hit cancel when finished');
    files = getFiles('*.txt','Choose a traces file:');
end

nFiles = numel(files);

if nFiles == 0,
    disp('No files specified, exiting.');
    return;
end



figure;

f = {};

% Load FRET histograms
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

    plot( fretaxis, pophist, 'Color',colors(i,:), 'LineWidth',2 );
    hold on;
    
    if isempty(f),
        f = fretaxis;
    end
    
    f = [ f pophist ];
    
%     area( fretaxis+i*0.003+(contour_bin_size/2), pophist, 'LineStyle', 'none', 'FaceColor', colors(i,:) );
%     alpha(0.2);
    
    % Setup visual style
end

hold off;
ylabel( 'Percent of total time' );

xlabel( 'FRET Efficiency' );
xlim( [0.1 1.0] );


if nargin>0,
    legend( titles );
end


save( 'pophist.txt', 'f', '-ASCII' );


end









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

