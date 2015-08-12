function plotDwelltimes( dwtfilename, timeConstants )
% PLOTDWELLTIMES
%
%   PLOTDWELLTIMES( DWT_filename, timeConstants )
%   plots the time spent in each FRET class as given in a dwell-time file.
%   timeConstants is a NxN cell array (N FRET classes), on for each 
%   type of possible transition.  Each element is a Mx2 array
%   (M aggregated states) with time constants in the first row
%   and amplitudes in the second row.  If an element is simply
%   a single value, it is taken to be the time-constant and the
%   amplitude is assumed to be 100%.  Elements on the diagonal
%   are ignored.

% A NOTE ON SMOOTH BINNING FROM QUB:
% the histogram is binned logarithmically, so some of the
% smallest bins contain no exact multiples of Δt. We simulated
% data and recorded each event's true duration as well as sampled
% duration, sampling at Δt. We observed that events with sampled
% duration kΔt have true duration distributed between
% (k − 1)Δt and (k + 1)Δt with a peak at kΔt and linear fall-off.
% This option performs that redistribution among neighboring bins. 



% if file not specified, prompt user
if nargin<1,
    [datafile,datapath] = uigetfile({'*.dwt'},'Choose a dwell-time file:');
    if datafile==0, return; end  %user hit "cancel"

    dwtfilename = [datapath datafile];
    disp(dwtfilename);
end

%% Old version

% Load the dwell times
[dwells,sampling,offsets,model] = loadDWT(dwtfilename);

nTraces  = numel(dwells);
nClasses = numel(model)/2;

% Build list of dwell times in each class
dwellc    = cell(nClasses,1);

for i=1:nTraces,
    
    classes = dwells{i}(:,1);
    times   = dwells{i}(:,2).*sampling;
    
    for j=1:nClasses,
        dwellc{j} = [dwellc{j} ; times(classes==j)];
    end
    
end

grey = repmat(0.75, 1,3);

% Build log-space histograms of dwell-times
% xaxis = power(10, start:step:end )
nbins = 40;
xaxis = logspace(1, 5, nbins);
xaxis = [unique( floor( xaxis/sampling )*sampling )];

% Create normalization factor for bin size (log-space)
lxa = log10(xaxis);
dlx = diff( [lxa(1:end-1) ; lxa(2:end)] );
dlx = [dlx dlx(end)];


% Build log-space histograms of dwell-times
% dwellhist = zeros(nClasses,nbins);
for i=1:nClasses,
    subplot(2,nClasses,i); cla;
    
    d = double( dwellc{i} );
    [N,NX] = hist( d, xaxis );
    N_norm = N./dlx;  %normalize by log-space bin size
    N_norm = N_norm./(sum(N_norm));  %normalize to 1
    
    h = bar( NX,N_norm, 'hist' );
    set(h, 'FaceColor',grey);
    set(gca, 'XScale', 'log');
    %xlabel( 'Dwell Time (log_{10} ms)' );
    ylabel( 'Events/total' );
    xlim( [30 10^5] );
    
    hold on;
    
end

% Plot kinetic fits (Sigworth & Sine, linear ordinate)
% for i=1:nClasses,
%     subplot( 1, nClasses, i );
%     
%     % g(x) is the Sigworth&Sine function, where x is
%     % in natural-log space.  We then plot it in log10 space.
%     tau = timeConstants(i);
%     x0 = log( tau );
%     x = 15 * exp( (0:40)*0.21 );
%     g = exp( log(x)-x0 - exp(log(x)-x0) );
%     g = g./sum(g);
% 
%     plot( x, g,'r-', 'LineWidth',3 );
% end


%% NEW VERSION = 
% Plot time in state before a transition to the other state occurs.

% Load the dwell times
[dwells,sampling,offsets,model] = loadDWT(dwtfilename);

nTraces  = numel(dwells);
nClasses = numel(model)/2;

% Build list of dwell times in each class
dwellc    = cell(nClasses,1);

for i=1:nTraces,
    
    classes = dwells{i}(:,1);
    times   = dwells{i}(:,2).*sampling;
    nDwells = numel(classes);
    
    start = find(classes~=1,1,'first'); % Skip initial dwells in dark state
    curClass = classes(start);  %class last observed
    curTime  = 0;               %total time so far spent in curState
    
    for j=start:nDwells
        
        % Ignore dwells in the dark states
        if classes(j) == 1,
            curTime = curTime + times(j);
        
        % If we are still in the same state, continue summing dwell times
        elseif classes(j) == curClass,
            curTime = curTime + times(j);
        
        % If a new state is encountered, add total time in the previous
        % state and add it to the distribution.
        else
            dwellc{curClass} = [ ...
                dwellc{curClass} ; curTime ];
            
            % Reset to new internal state
            curClass = classes(j);
            curTime  = times(j);
        end
        
    end %for each dwell
        
end %for each trace

for i=1:nClasses,
    disp( numel(dwellc{i}) );
end

grey = repmat(0.75, 1,3);

% Build log-space histograms of dwell-times
% xaxis = power(10, start:step:end )
nbins = 40;
xaxis = logspace(1, 5, nbins);
xaxis = unique( floor( xaxis/sampling )*sampling );

% Create normalization factor for bin size (log-space)
lxa = log10(xaxis);
dlx = diff( [lxa(1:end-1) ; lxa(2:end)] );
dlx = [dlx dlx(end)];


% Build log-space histograms of dwell-times
% dwellhist = zeros(nClasses,nbins);
for i=2:nClasses,
    subplot(2,nClasses,nClasses+i); cla;
    
    d = double( dwellc{i} );
    [N,NX] = hist( d, xaxis );
    N_norm = N./dlx;  %normalize by log-space bin size
    N_norm = N_norm./(sum(N_norm));  %normalize to 1
    
    h = bar( NX,N_norm, 'hist' );
    set(h, 'FaceColor',grey);
    set(gca, 'XScale', 'log');
    xlabel( 'Dwell Time (log_{10} ms)' );
    ylabel( 'Events/total' );
    xlim( [30 10^5] );
    
    hold on;
    
end


% Plot kinetic fits (Sigworth & Sine, linear ordinate)
% for i=1:nClasses,
%     subplot( 1, nClasses, i );
%     
%     % g(x) is the Sigworth&Sine function, where x is
%     % in natural-log space.  We then plot it in log10 space.
%     tau = timeConstants(i);
%     x0 = log( tau );
%     x = 15 * exp( (0:40)*0.21 );
%     g = exp( log(x)-x0 - exp(log(x)-x0) );
%     g = g./sum(g);
% 
%     plot( x, g,'r-', 'LineWidth',3 );
% end






end %function plotDwelltimes
