function frethist = makecplot( data, options)
% MAKECPLOT   creates _hist.txt FRET histogram file

% Load data
if ischar(data)
    data_filename = data;
    [d,a,fret] = loadTraces( data_filename );
else
    fret = data;
end

% Load options
if nargin<2,
    options.contour_bin_size = 0.03;
end

if ~isfield(options,'options.pophist_offset'),
    options.options.pophist_offset = 0;
end


% Cut off first few frames to get rid of gain drift.
fret = fret( :, (1+options.pophist_offset:end) );
[Nmol len] = size(fret);

% Axes for histogram includes all possible data
time_axis = 1:len;
fret_axis = -0.1:options.contour_bin_size:1.2;

% Initialize histogram array, setting the time step in the first row,
% and the FRET bins in the first column. This is done for import into
% Origin.
frethist = zeros( length(fret_axis)+1, length(time_axis)+1 );
frethist(1,2:end) = time_axis;
frethist(2:end,1) = fret_axis';

frethist(2:end,2:end) = hist( fret, fret_axis  );


% Save plots to file
if exist('data_filename','var'),
    histfile=strrep(data_filename,'.txt','_hist.txt');
    dlmwrite(histfile,frethist,' ');

    % Also save a normalized histogram.
    frethistn = frethist;
    frethistn(2:end,2:end) = frethistn(2:end,2:end)/Nmol;
    histfile=strrep(data_filename,'.txt','_normhist.txt');
    dlmwrite(histfile,frethistn,' ');
end


end %function makecplot
