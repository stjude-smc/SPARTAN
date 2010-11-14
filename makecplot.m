function frethist = makecplot( data, options)
% MAKECPLOT   creates _hist.txt FRET histogram file

% Load data
if ischar(data)
    data_filename = data;
    [d,a,fret] = loadTraces( data_filename );
    
    if size(d,1)<1,
        error('File is empty: %s',data_filename);
    end
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

% Remove traces with NaN values
bad = isnan( fret(:) );
if any(bad),
    fret(bad) = 0;
    warning('NaN values found in FRET data. Converting to zeros.');
end

% Axes for histogram includes all possible data
time_axis = 1:len;
fret_axis = -0.1:options.contour_bin_size:1.2;

% Initialize histogram array, setting the time step in the first row,
% and the FRET bins in the first column. This is done for import into
% Origin.
frethist = zeros( length(fret_axis)+1, length(time_axis)+1 );
frethist(1,2:end) = time_axis;
frethist(2:end,1) = fret_axis';

% WARNING: if the fret variable is too large (>6M datapoints?) hist will fail
% and return all zeros. As a hack, the traces are split into smaller batches
% when the total size is large (>2M datapoints). FIXME
% if numel(fret)>2e6,
    nFrames=size(fret,2);
    for i=1:1000:nFrames,
        t_range = i:min(i+1000-1,nFrames);
        frethist(2:end,1+t_range) = hist( fret(:,t_range), fret_axis  );
    end
% else
%     frethist(2:end,2:end) = hist( fret, fret_axis  );
% end


% Save plots to file
if exist('data_filename','var'),
    [p,n,e]  = fileparts(data_filename);
    histfile = [p filesep n '_hist.txt'];
    dlmwrite(histfile,frethist,' ');

    % Also save a normalized histogram.
    frethistn = frethist;
    frethistn(2:end,2:end) = frethistn(2:end,2:end)/Nmol;
    histfile = [p filesep n '_normhist.txt'];
    dlmwrite(histfile,frethistn,' ');
end


end %function makecplot
