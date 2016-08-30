function varargout = fluorcorr(filenames,minFret)
% FLUORCORR:  Cross correlation of donor and acceptor fluorescence traces.
% 
%   CORR = FLUORCORR( FILENAMES )  loads fluorescence traces from the specified
%   filenames (a cell array or a string for a single file), calculates the
%   cross-correlation of donor and acceptor fluorescence traces, and takes an
%   average across all traces in the file (is that correct?). Only the region
%   before acceptor or donor photobleaching is considered so that any
%   anticorrelation between donor and acceptor traces is due to dynamics and not
%   photophysics. Blinking may still contribute (FIXME).
%
%   CORR is an array of files in columns and lag times in rows. The first column
%   contains the X-axis (lag times in frames). The output is also saved as a
%   file for easy importing into Origin called "fluorcorr.txt".
%
% See also: traceStat.
%

% TODO: consider removing some outlier traces that have many blinking events,
% since we are only removing a fraction of them with the threshold. Also
% consider dilating the binary mask that defines when the dyes are alive to
% ensure the rising and falling edges do not contribute to xcorr.
% Should test all this with a simulation with very fast, but visible, blinking.
%


%% CONSTANTS
maxlag = 50;   %maxmim lag interval over which to calculate cross correlation.
                %large lags have fewer datapoints, making them poorly defined (noisy).

% Threshold to remove acceptor (and donor) dark states. Set this value as high
% as possible but below all real (conformational) FRET states. Blinking can
% contribute significantly to xcorr, so it must be removed for analysis.
if nargin<2 || isempty(minFret),
%     minFret = 0.15;  %default
    minFret = 0.18;  %trying to remove at least some of the blinking with LeuT.
end

% Stolen from frethistComparison. Maybe it should be in cascadeConstants.
% FIXME: like frethistComparison, have a fallback for a large number of files.
colors = [ 0      0      0    ; ...  % black
           0.75   0      0.75 ; ...  % purple
           0      0.75   0.75 ; ...  % cyan
           0      0.5    0    ; ...  % green
           0.75   0.75   0    ; ...  % yellow
           1      0      0    ; ...  % red
           0.6    0      0    ; ...     % dark red
           0.6    0      0    ; ...     % dark red
           0.6    0      0    ; ...     % dark red
           0.6    0      0    ];     % dark red
               
%% Load the data
if nargin<1,
    filenames = getFiles;
end
if ~iscell(filenames),
    filenames = {filenames};
end
nFiles = numel(filenames);

if nFiles==0, return; end

figure; hold on;
if nFiles > 7,
    colors = zeros(nFiles,3);
    interval = (1/(nFiles-1));
    colors(:,1) = 0:interval:1;
    colors(:,3) = 1:-interval:0;
end
set(gca,'ColorOrder',colors);


%% Calculate cross-correlation on each file.
output = zeros( numel(filenames), maxlag+1 );

for f=1:numel(filenames)
    % Load the trace data from file.
    data = loadTraces(filenames{f});
    data.truncate(100);
    nTraces = 0;  %number of traces actually used for the calculation.
    
    % For each trace, calculate xcorr b/t donor and acceptor.
    for i=1:data.nTraces,
        % Use only the region prior to photobleaching so only dynamics between
        % non-zero FRET states is considered.
        % NOTE: xcorr time lag may be distorted if this filter removes a 
        % significant number of datapoints because it is skipping ahead of time!
        idxAlive = data.fret(i,:)>minFret;
        
        % Remove the rising and falling edges of each of these blinking events.
        % Apparently this has almost no effect.
        idxAlive = imerode(idxAlive, ones(1,3));
        
        % Calculate cross correlation for traces that are long enough to fill
        % the maximum lag, and then some.
        if sum(idxAlive) >= maxlag*2,
            nTraces = nTraces+1;
            D = data.donor(i,idxAlive);
            A = data.acceptor(i,idxAlive);
            
            [C,lags] = xcorr( D-mean(D), A-mean(A), maxlag, 'coeff' );            
            output(f,:) = output(f,:) + C( find(lags==0):end );
        end
    end
    output(f,:) = output(f,:)/nTraces;
    
    % Add the average cross-correlation curve of this file to the plot.
    plot( (0:maxlag)', output(f,:), 'Color',colors(f,:), 'LineWidth',3 );
end


xlabel('Lag time (frames)');
ylabel('Cross-correlation');
legend( trimtitles(filenames) );
xlim([0 maxlag]);



%% Save the result for easy import into Origin.
output = [ (0:maxlag)' output' ];
save('fluorcorr.txt','output','-ASCII');

if nargout>0,
    varargout{1} = output;
end

end




