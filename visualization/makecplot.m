function [hist2d,units] = makecplot( input, options )
% MAKECPLOT   Creates a contour plot of FRET values over time.
%
%   [HIST] = MAKECPLOT(INPUT, OPTIONS) creates a FRET-time contour
%   histogram (HIST) from the INPUT Traces object or .traces filename.
%   HIST: first row is time (seconds), first column is FRET value of each bin.
%   OPTIONS is a struct array for additional options -- see cascadeConstants.m.
%
%   See also: makeplots, cplot.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


% Load traces and select appropriate FRET field.
if ischar(input),
    data = loadTraces(input);
elseif isa(input,'TracesFret') || isstruct(input)
    data = input;
else
    error('Invalid trace data input');
end

% Select FRET data and remove initial frames if requested.
fret = data.(options.fretField)(:,1+options.pophist_offset:end);
if options.frameAxis
    time_axis = options.pophist_offset + (1:data.nFrames);
    units = 'frames';
else
    time_axis = data.time(1+options.pophist_offset:end)/1000;
    units = 's';
end

% Remove traces with NaN values
bad = isnan( fret(:) );
if any(bad),
    fret(bad) = 0;
    warning('NaN values found in FRET data. Converting to zeros.');
end

% Initialize histogram array, setting the time step in the first row,
% and the FRET bins in the first column. This is done for import into
% Origin.
fret_axis = options.fret_axis;
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

% Normalize. "to_max" option normalizes to max among several plots.
if isfield(options,'cplot_normalize_to_max') && options.cplot_normalize_to_max,
    hist2d(2:end,2:end) = hist2d(2:end,2:end)/options.cplot_normalize_to_max;
else
    hist2d(2:end,2:end) = hist2d(2:end,2:end)/size(fret,1);
end

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


end %function makecplot
