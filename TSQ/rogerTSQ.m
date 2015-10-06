function rogerTSQ(fnames)
%rogerTSQ  Fluorophore performance statistics (single-color)
%
%   rogerTSQ(FILES) calculates trace statistics relevant for the evaluation
%   of the photophysical behavior of the dyes in each filename in the
%   cell array FILES. These will be printed in the command window and
%   displayed in a summary plot. The analysis procedure uses automatic
%   hidden Markov modeling and idealization of fluorescence traces.
%
%   Intensity: mean total fluorescence intensity before photobleaching.
%   SNR:       signal-to-backround noise ratio
%   Ton:       exponential lifetime of the fluorescent state between blinks.
%   Toff:      exponential lifetime in the non-fluorescent (blinking) state
%   Total Ton: total time in the fluorescent state before photobleaching.
%   Yield:     Total photon yield (approximately Intensity*Total Ton).
%   pt:        percentage of total time in the fluorescent state.
%
%   NOTE: photobleaching is marked as the last large drop in total
%   fluorescence intensity. Traces that are not photobleached at the end of
%   the movie may be excluded from analysis.

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% TODO: add errors bars for stats, use stretched exponential fit for dwell times
% to account for heterogeneity. The stretch factor could be plotted also as an
% error to show the spread in the behavior. Beta closer to 0 would give high
% errors, close to 1 would give no error. (1-Beta)*tau

% TODO: ask what fluorescence channel(s) to analyze.



%% Options

% Save traces rejected from analysis in a separate file for testing.
SAVE_REJECTED = false;


%%
%-------------------------------------------------------------
% 1) Request filenames from user.
if nargin<1,
    fnames = getFiles('*.rawtraces');
    if isempty(fnames), return; end
end
nFiles = numel(fnames);

% Define Markov model for blinking events.
% FIXME: should we just create a model?
model = qub_loadModel;
if isempty(model), return; end  %user hit cancel.

model.fixSigma = [1 0];
skmParams.seperately = 1;
skmParams.quiet = 1;

mu    = model.mu';
sigma = model.sigma';
onState  = find( mu==max(mu) );
offState = find( mu==min(mu) );


% prep output variables
intensity    = zeros( nFiles,1 );
intensityStd = zeros( nFiles,1 );
SNRs         = zeros( nFiles,1 );
SNRsStd      = zeros( nFiles,1 );
Ton          = zeros( nFiles,1 );
Toff         = zeros( nFiles,1 );
totalTon     = zeros( nFiles,1 );
totalTonStd  = zeros( nFiles,1 );
yield        = zeros( nFiles,1 );
yieldstd     = zeros( nFiles,1 );

names        = cell( nFiles,1 );  %plain-text name of each file
dwtFilename  = cell( nFiles,1 );


hf = figure;
set(hf,'Units','normalized');
set(hf,'Position', [0.16  0.29 0.67 0.4] );


%%
for i=1:nFiles,
    
    %-------------------------------------------------------------
    % 1) Load data from each file and calculate stats. 
    data = loadTraces( fnames{i} );
    sampling = data.sampling/1000;
    
    [p,f] = fileparts( fnames{i} );
    names{i} = strrep(f,'_',' ');
    basename = fullfile(p,f);
    
    
    %-------------------------------------------------------------
    % 2) Pick traces according to defined criteria if not already filtered.
    stats = traceStat( data );
    
    % Fit the background distribution to find a good cutoff to remove
    % traces with low-intensity resurrection (not well-handled by HMM).
    % NOTE: this may exclude traces that blink but do not photobleach.
    bins = 0:70;
    bg = [stats.bg];
    bg = bg( bg>bins(1) & bg<bins(end) );
    N = hist( bg, bins );
    f = fit( bins', N', 'gauss1', 'StartPoint',[15,median(bg),0.5*median(bg)] );
    cutoff = f.b1 + 4*f.c1; %4 deviations from mean.

    % Define selection criteria
    criteria.min_snr = 10;  %FIXME: this seems too high??
    criteria.eq_overlap = 0;
    criteria.max_bg = cutoff;

    % Select traces according to criteria defined above.
    indexes = pickTraces( stats, criteria );
    data.subset(indexes);
    stats = stats(indexes);
    
    
    %-------------------------------------------------------------
    % 3) Calculate signal statistics
    
    % Get total intensity distribution parameters.
    % Histogram fitting is best with clear, symmetric distributions.
    t = [stats.t];
    [histdata,bins] = hist( t(t>0), 40 );
    histdata = 100*histdata/sum(histdata);  %normalize
%     f = fit( bins',histdata', 'gauss1' );
%     intensity(i) = f.b1;
%     intensityStd(i) = sqrt(f.c1/2);
    intensity(i) = median(t);
    intensityStd(i) = std( bootstrp(1000,@median,t) );

    subplot(2,5,1); hold on;
    plot( bins, histdata );
    xlabel('Intensity (photons)');  ylabel('Counts (%)');
    
    
    snr = [stats.snr_s];
    [histdata,bins] = hist( snr(snr>0), 40 );
    histdata = 100*histdata/sum(histdata);  %normalize
%     f = fit( bins',histdata', 'gauss1' );
%     SNRs(i) = f.b1; %mean
%     SNRs(i) = sqrt(f.c1/2); %std
    SNRs(i) = median(snr);
    SNRsStd(i) = std( bootstrp(1000,@median,snr) );

    subplot(2,5,2); hold on;
    plot( bins, histdata );
    xlabel('SNR');  ylabel('Counts (%)');
    
    
    %-------------------------------------------------------------
    % 4) Idealize the intensity data to a two-state model using SKM. 
    % Total intensity is scaled to 1 on average, analogous to cy5forQuB.
    dwt = skm( data.total/intensity(i), data.sampling, model, skmParams );
    
    % FIXME: consider removing traces with no ON-state dwells.
    
    %-------------------------------------------------------------
    % 5) Remove traces that are poorly idealized or are significant outliers.
    % A) Remove any trace with >2 stdev transitions/total lifetime.
    %donorlife = [stats.donorlife];
    nTransitions = cellfun('size',dwt,1)-1;
    %donorlife = reshape( donorlife, numel(donorlife),1 );
    nTransitions = reshape( nTransitions, numel(nTransitions),1 );
    
    transRate = nTransitions;%./donorlife;
    meanTrans = mean(transRate);
    stdTrans  = std(transRate);
    
    if stdTrans>0,
        selected = transRate < (meanTrans + 2*stdTrans);
        nRejected = data.nTraces - sum(selected);
    else
        selected = true(size(transRate));
        nRejected = 0;
    end
    
    fprintf('File %d: Removed %d traces (%.1f%%)\n', i, nRejected, 100*nRejected/data.nTraces);
  
    
    %-------------------------------------------------------------
    % 6) Save idealization.
    
    % Save selected traces and idealization.
    saveTraces( [basename '_auto.traces'], data.getSubset(selected) );
    
    dwtFilename{i} = [basename '.qub.dwt'];
    offsets = data.nFrames*((1:sum(selected))-1);
    saveDWT( dwtFilename{i}, dwt(selected), offsets, [mu sigma], data.sampling );
    
    % Save rejected traces and idealization.
    if SAVE_REJECTED,
        saveTraces( [basename '_rejected.traces'], data.getSubset(~selected) );

        offsets2 = data.nFrames*((1:sum(~selected))-1);
        saveDWT( [basename '_rejected.dwt'], dwt(~selected), offsets2, ...
                                                        [mu sigma], data.sampling );
    end
    
         
    %-------------------------------------------------------------   
    % 6) Calculate total lifetimes      
    onTimes  = [];
    offTimes = [];
    
    % Collect list of dwell-times in each class  
    for j=1:numel(dwt),
        classes  = dwt{j}(:,1);
        times    = dwt{j}(:,2).*sampling;
        
        % Remove final off-state dwell, which may be photobleached state.
        if classes(end)==offState,
            classes = classes(1:end-1);
            times   = times(1:end-1);
        end
        
        onTimes  = [onTimes  ; times(classes==onState ) ];
        offTimes = [offTimes ; times(classes==offState) ];
    end
    
    % Total time in the ON state for each trace.
    dwt = dwt(selected);
    idl = dwtToIdl( dwt, data.nFrames, offsets, data.nTraces );
    totalOn = sum(idl==onState,2).*sampling;
    
    % Calculate average total time on/off
    Ton(i)  = mean( onTimes  );
    Toff(i) = mean( offTimes );
    totalTon(i)  = mean( totalOn );
    totalTonStd(i) = std( bootstrp(1000,@mean,totalOn) );
    
    % Calculate photon yield
    photons = sum( data.total.*(idl==onState), 2 );
    yield(i)    = mean( photons );
    yieldstd(i) = std( bootstrp(1000,@mean,photons) );
    
    
    %-------------------------------------------------------------
    % 7) Display state lifetime histograms.    
    dwellaxis = (0:1:data.nFrames)*sampling;  %fixme
    
    % ON times
    subplot(2,5,3); hold on;
    survival(onTimes,dwellaxis);
    xlabel('Time ON (s)');  ylabel('Counts (%)');
    
    % OFF times
    subplot(2,5,4); hold on;
    survival(offTimes,dwellaxis);
    xlabel('Time OFF (s)');  ylabel('Counts (%)');
    
    % Total time ON
    subplot(2,5,5); hold on;
    survival(totalOn,dwellaxis);
    xlabel('Total time ON (s)');  ylabel('Counts (%)');
end


% Ask the user to give the files specific names.
% prompts = cellfun( @(n)sprintf('File %d',n), num2cell(nFiles), 'UniformOutput',false );
% answer = inputdlg( prompts, 'Enter short file titles',1, names );
% if ~isempty(answer)
%     names = answer;
% end

subplot(2,5,5); legend(names);
% FIXME: move into the open space?


%% Make bar graphs to summarize the results

subplot(2,5,5+1); hold on;
bar( 1:nFiles, intensity, 'r' );
errorbar( 1:nFiles, intensity, intensityStd, '.k' );
ylabel('Intensity (photons)');
xlim([0.35 nFiles+0.65]);

subplot(2,5,5+2); hold on;
bar( 1:nFiles, SNRs, 'r' );
errorbar( 1:nFiles, SNRs, SNRsStd, '.k' );
ylabel('Signal-noise ratio');
xlim([0.35 nFiles+0.65]);

subplot(2,5,5+3);
bar( Ton, 'r' );
ylabel('Time ON (s)');
xlim([0.35 nFiles+0.65]);

subplot(2,5,5+4); hold on;
bar( totalTon, 'r' );
errorbar( 1:nFiles, totalTon, totalTonStd, '.k' );
ylabel('Total time ON (s)');
xlim([0.35 nFiles+0.65]);

subplot(2,5,5+5); hold on;
bar( 1:nFiles, yield, 'r' );
errorbar( 1:nFiles, yield, yieldstd, '.k' );
ylabel('Total photon yield');
xlim([0.35 nFiles+0.65]);



%%

% Calculate dwell times by fitting to exponential decays..
% FIXME
% params.useCorrectedDwelltimes = 0; %don't merge blinks)
% params.fitSingle = true;
% params.quiet = true;
% Ton = lifetime_exp( dwtFilename, params );
% Ton = Ton(:,2);

%pt = percentTime(dwtFilename);
%pt = pt(:,2);

% Save results to a file that can be plotted as bar graphs in Origin.
% FIXME: use a Matlab function for this...
[f,p] = uiputfile('tsqStats.txt','Choose a filename to save the results.');
if f,
    fid = fopen( fullfile(p,f), 'w' );

    fwrite(fid, sprintf('Name\tIntensity (photons)\tIntensity stdev\tSNR\tSNR stdev\ttON (sec)\tTotal tON (sec)\r\n') );

    for i=1:nFiles,
        fwrite( fid, sprintf('%s\t%.0f\t%.0f\t%.1f\t%.1f\t%.2f\t%.2f\r\n', ...
                names{i}, intensity(i),intensityStd(i), SNRs(i),SNRsStd(i), Ton(i), totalTon(i) )   );
    end
    
    fclose(fid);
end



end %function rogerTSQ



function histdata = survival(times,dwellaxis)

histdata = histc( times, dwellaxis );
histdata = sum(histdata) - cumsum(histdata);  %survival plot
histdata = histdata/histdata(1);  %normalize

plot( dwellaxis, histdata );
    
end



