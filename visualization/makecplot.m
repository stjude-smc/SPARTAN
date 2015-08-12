function [hist2d] = makecplot( input, options )
% MAKECPLOT   Creates a contour plot of FRET values over time.
%
%   [HIST] = MAKECPLOT( FRET )
%   Sums FRET values from FRET (traces in rows) into a histogram at each
%   point in time (HIST).
%
%   [HIST] = MAKECPLOT( FILENAME )
%   Sums FRET values from data loaded from FILENAME. The histogram data is
%   also saved in a file with the extension "_normhist.txt".
%
%   [...] = MAKECPLOT( ..., options )
%   A structure containing options can be given to specify how the
%   histograms should be made. These are listed below in makeplots.m.
%   FIXME: These actually need some work for consistency and setting all
%   the defaults in cascadeConstants.m
%

% Load data
if ischar(input)
    data_filename = input;
    data = loadTraces( data_filename );
    fret = data.fret;
    
    if size(data.donor,1)<1,
        error('File is empty: %s',data_filename);
    end
elseif isstruct(input)
    fret = input.fret; %assuming traces data structure
elseif isnumeric(input),
    fret = input;
else
    error('Unknown input data parameter type');
end

% Load options
if nargin<2,
    options = constants.defaultMakeplotsOptions;
end

if ~isfield(options,'pophist_offset'),
    options.pophist_offset = 0;
end


% Cut off first few frames to get rid of gain drift.
fret = fret( :, 1+options.pophist_offset:end );
[Nmol,len] = size(fret);

% Remove traces with NaN values
bad = isnan( fret(:) );
if any(bad),
    fret(bad) = 0;
    warning('NaN values found in FRET data. Converting to zeros.');
end

% Axes for histogram includes all available data.
time_axis = 1:len;
fret_axis = options.fret_axis;

% Initialize histogram array, setting the time step in the first row,
% and the FRET bins in the first column. This is done for import into
% Origin.
hist2d = zeros( length(fret_axis)+1, length(time_axis)+1 );
hist2d(1,2:end) = time_axis;
hist2d(2:end,1) = fret_axis';

% Calculate histograms over each time bin. NOTE: if the fret variable is
% too large (>6M datapoints?) hist will fail and return all zeros.
% As a hack, the traces are split into smaller batches.
nFrames = size(fret,2);

for i=1:1000:nFrames,
    t_range = i:min(i+1000-1,nFrames);
    
    % Pad the data with NaN values (which do not contribute to the histogram) so
    % that hist() acts the same if there is only one trace. Without this, hist
    % with produce a single histogram for the entire trace.
    padded = [  fret(:,t_range) ;  NaN(1,numel(t_range))  ];
    
    hist2d(2:end,1+t_range) = hist( padded, fret_axis  );
end


% Normalize the histogram to the total number of traces. This way all plots
% can use the same scale.
hist2d(2:end,2:end) = hist2d(2:end,2:end)/Nmol;


% Optional: remove traces as they photobleach (drop to zero-FRET). This way the
% plot looks the same as traces bleach -- there are simply fewer of them
% contributing.
if isfield(options,'cplot_remove_bleached') && options.cplot_remove_bleached,
    % Set bleached bins to zero.
    % FIXME: this threshold should be defined in cascadeConstants.
    h = hist2d(2:end,2:end);
    h( fret_axis<=0.15, : ) = 0;

    % Renormalize so all columns sum to unity again.
    hist2d(2:end,2:end) = bsxfun( @rdivide, h, sum(h) );
end


% Save the histogram to file. Only do this if user doesn't request output
% arguments. What else could they want?
if nargout==0 && exist(data_filename,'var'),
    [p,n] = fileparts(data_filename);
    hist_fname = fullfile(p, [n '_hist.txt']);
    dlmwrite(hist_fname,hist2d,' ');
end


end %function makecplot
