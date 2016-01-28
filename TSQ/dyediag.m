function [output,errors,hfig] = dyediag(fnames, modelInput)
%dyediag  Fluorophore performance statistics (single-color)
%
%   dyediag(FILES) calculates trace statistics relevant for the evaluation
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
%   If multiple files are selected in a directory, they are assumed to be
%   from the same condition and are analyzed together. In that case, error
%   bars are the standard deviation across the movies. If only one file is
%   selected for a condition, error bars are estimated from bootstrap
%   samples.
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


narginchk(0,2);
nargoutchk(0,3);


%% Options

% Save traces rejected from analysis in a separate file for testing.
SAVE_REJECTED = false;


%%
%-------------------------------------------------------------
% 1) Request filenames from user.
% fnames is a cell array, each element of which is a cell array of
% filenames that are presumably from the same condition. This requires the
% user to ensure all files from a particular conidtion are in the same
% folder, as would be expected if everything is from one day.
if nargin<1,
    fnames = getFileGroups('*.rawtraces');
end
nFiles = numel(fnames);  %actually the number of file groups/conditions.
if isempty(fnames), return; end


% Load model for detecting blinks. States: 1) bleached, 2) blinking, 3) ON.
if nargin>=2,
    if isa(modelInput,'QubModel')
        model = modelInput;
    elseif isstruct(modelInput) || ischar(modelInput),
        model = QubModel(modelInput);
    else
        error('Invalid model input. Must be a QubModel object or file name');
    end
else
    % Define a standard internal model
    model.p0    = [0.01 0.01 0.98]';
    model.class = [1 1 2];
    model.mu    = [0 1];
    model.sigma = [0.1 0.1]; %fixme?
    model.rates = [0        0       0
                   0.01     0       1    %bleaching, ressurection rates.
                   0        1       0];  %bleaching, blinking rate (s-1)
    model = QubModel(model);
end

model.fixSigma = [1 0];
skmParams.seperately = 1;
skmParams.quiet = 1;

mu    = to_col(model.mu);
sigma = to_col(model.sigma);
onState  = find( mu==max(mu) );
offState = find( mu==min(mu) );


% Prep output variable list. The order here also determines the order of 
% the columns in the output file.
z = { zeros(1,nFiles) };
output = struct('intensity',z, 'SNRs',z, 'Ton',z, 'Toff',z, 'totalTon',z, 'yield',z);
errors = output;

% Plain-text column names for axes labels and file output.
colnames = {'Intensity (photons)','SNR', 'Time ON (s)', 'Time OFF (s)', ...
            'Total time ON (s)','Photon yield'};
             
names = cell( nFiles,1 );  %plain-text name of each file


constants = cascadeConstants;
hfig = figure( 'Name',sprintf('Dye Diagnostics (ver. %s)', constants.version), ...
               'Units','normalized', 'Position', [0.16  0.29 0.67 0.4] );
ax = zeros(2,5);  %subplot axes


%%
for i=1:nFiles,
    
    %-------------------------------------------------------------
    % 1) Load data from each file and combine into one large dataset.
    condition_files = fnames{i};
    condition_data  = cell(  numel(condition_files), 1 );
    condition_idx   = [];
    
    [p,f] = fileparts( condition_files{1} );
    names{i} = strrep(f,'_',' ');
    basename = fullfile(p,f);
    
    for j=1:numel(condition_files),
        condition_data{j} = loadTraces( condition_files{j} );
    end
    
    % Combine data from all files into one larger Traces object.
    if numel(condition_data)>1,
        data = combine( condition_data{:} );
    else
        data = condition_data{1};
        fprintf('Only one file for condition #%d. Using bootstrap sampling for error bars\n',i);
    end
    
    % Keep track of which file each trace came from.
    for j=1:numel(condition_data),
        condition_idx = [condition_idx ; repmat(j,condition_data{j}.nTraces,1)];
    end
    clear condition_data;
    
    sampling = data.sampling/1000;
    
    
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

    % Remove any extreme traces (1%) so histograms are reasonably scaled.
    t   = sort( [stats.t]     );
    snr = sort( [stats.snr_s] );
    tMax   = t( floor(0.99*data.nTraces)-1 );
    snrMax = snr( floor(0.99*data.nTraces)-1 );
    
    % Define selection criteria
    criteria.min_snr = 10;
    criteria.eq_overlap = 0;
    criteria.max_bg = cutoff;
    criteria.min_t = 0;
    criteria.max_t = tMax;
    criteria.min_snr_s  = 0;
    criteria.max_snr_s = snrMax;

    % Select traces according to criteria defined above.
    [selected,stats] = pickTraces( stats, criteria );
    if numel(selected)<1, continue; end
    
    data.subset(selected);
    condition_idx = condition_idx(selected);
    
    
    %-------------------------------------------------------------
    % 3) Calculate signal statistics
    
    % Total intensity distributions
    t = [stats.t];
    output.intensity(i) = median(t);
    errors.intensity(i) = stdbyfile(t,condition_idx,@median);
    
    [histdata,bins] = hist( t, 40 );
    histdata = 100*histdata/sum(histdata);  %normalize
    ax(1,1) = subplot(2,5,1); hold all;
    plot( bins, histdata );
    
    
    % Signal-to-noise over signal distributions
    snr = [stats.snr_s];
    output.SNRs(i) = median(snr);
    errors.SNRs(i) = stdbyfile(snr,condition_idx,@median);
    
    [histdata,bins] = hist( snr, 40 );
    histdata = 100*histdata/sum(histdata);  %normalize
    ax(1,2) = subplot(2,5,2); hold all;
    plot( bins, histdata );
    
    
    %-------------------------------------------------------------
    % 4) Idealize the intensity data to a two-state model using SKM. 
    % Total intensity is scaled to 1 on average, analogous to cy5forQuB.
    % Note: model re-estimation includes the photobleached state, which may
    % not be ideal for getting blinking kinetics.
    dwt = skm( data.total/output.intensity(i), data.sampling, model, skmParams );
    assert( numel(dwt)==data.nTraces, 'Idealization size mismatch' );
    
    % FIXME: consider a filter for brief events after photobleaching.
    % FIXME: consider removing traces with no ON-state dwells.
    
    
    %-------------------------------------------------------------
    % 5) Remove traces that are poorly idealized or are significant outliers.
    
    % Remove any trace with >2 stdev transitions/total lifetime.
    nTransitions = cellfun('size',dwt,1)-1;
    nTransitions = reshape( nTransitions, numel(nTransitions),1 );
    
    transRate = nTransitions;  %./[stats.donorlife];  %make it a rate
    meanTrans = mean(transRate);
    stdTrans  = std(transRate);
    
    if stdTrans>0,
        selected = transRate < (meanTrans + 2*stdTrans);
    else
        selected = true(size(transRate));
    end
    
    % Remove empty traces
    selected = selected & nTransitions>0;    
    nRejected = data.nTraces - sum(selected);
    
    fprintf('File %d: Removed %d traces (%.1f%%)\n', i, nRejected, 100*nRejected/data.nTraces);
  
    
    %-------------------------------------------------------------
    % 6) Save idealization.
    
    % Save rejected traces and idealization.
    if SAVE_REJECTED,
        saveTraces( [basename '_rejected.traces'], data.getSubset(~selected) );

        offsets2 = data.nFrames*((1:sum(~selected))-1);
        saveDWT( [basename '_rejected.dwt'], dwt(~selected), offsets2, ...
                                               [mu sigma], data.sampling );
    end
    
    % Save selected traces and idealization.
    data.subset(selected);
    dwt = dwt(selected);
    condition_idx = condition_idx(selected);    
    
    saveTraces( [basename '_auto.traces'], data );
    
    offsets = data.nFrames*( (1:data.nTraces)-1 );
    saveDWT( [basename '_auto.qub.dwt'], dwt, offsets, [mu sigma], data.sampling );
        
    
    %-------------------------------------------------------------   
    % 6) Calculate total lifetimes
    dwellaxis = (0:1:data.nFrames)*sampling;
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
    
    % Get total time on for each trace.
    idl = dwtToIdl( dwt, offsets, data.nFrames, data.nTraces );
    totalOn = sum(idl==onState,2).*sampling;
    output.totalTon(i) = mean(totalOn);   %expfit(totalOn, dwellaxis);
    errors.totalTon(i) = stdbyfile(totalOn,condition_idx,@mean);
    
    % Calculate average total time on/off
    output.Ton(i)  = mean( onTimes );
    %errors.Ton(i) = stdbyfile(onTimes,condition_idx,@mean);
    
    output.Toff(i) = mean( offTimes );
    %errors.Toff(i) = stdbyfile(offTimes,condition_idx,@mean);
    
    
    % Calculate photon yield
    photons = sum( data.total.*(idl==onState), 2 );
    output.yield(i) = mean( photons );
    errors.yield(i) = stdbyfile(photons,condition_idx,@mean);
    
    
    %-------------------------------------------------------------
    % 7) Display state lifetime histograms.
    ax(1,3) = subplot(2,5,3); hold all;
    survival(onTimes,dwellaxis);
    
    ax(1,4) = subplot(2,5,4); hold all;
    survival(offTimes,dwellaxis);
    
    ax(1,5) = subplot(2,5,5); hold all;
    survival(totalOn,dwellaxis);
    
    drawnow;
end

% Plot formatting
for i=1:size(ax,2),
    ylabel( ax(1,i), 'Counts (%)' );
    xlabel( ax(1,i), colnames{i} );
end

legend(names);


%% Make bar graphs to summarize the results

fieldIdx = [1:3,5:6];  %which fields to show
fields = fieldnames(output);

for i=1:numel(fieldIdx),
    fid = fieldIdx(i);
    statName = fields{fid};
    
    ax(2,i) = subplot(2,5,5+i); hold all;
    bar( 1:nFiles, output.(statName), 'r' );
    errorbar( 1:nFiles, output.(statName), errors.(statName)/2, '.k' );
    
    ylabel(colnames{fid});
    xlim([0.35 nFiles+0.65]);
end

xlabel(ax(2,1), 'Condition');
set(ax(2,:),'xtick',1:nFiles);
set(ax,'box','on');


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
[f,p] = uiputfile('tsqStats.txt','Choose a filename to save the results.');
if f,
    % Write header line with field labels and units (see top of file).
    fid = fopen( fullfile(p,f), 'w' );
    header = sprintf('%s\t','Name',colnames{:});  %FIXME add error headers
    fprintf(fid, '%s\n', header(1:end-1) );
    
    % Write data lines, one per file.
    fields = fieldnames(output);
    format = ['%s' repmat('\t%0.3f',[1 2*numel(fields)]) '\n'];

    for i=1:nFiles,
        line = cellfun( @(f)[output.(f)(i) errors.(f)(i)], fields, 'UniformOutput',false );
        fprintf( fid, format, names{i}, line{:} );
    end
    
    fclose(fid);
end



end %function dyediag



function stdev = stdbyfile(stat,idx,fcn)
% Select out all statistic values (stat) from each file as specified in
% idx, apply the fcn to get a mean (or median, etc) value, and then take
% the standard deviation across all such files.
% If only one file is given, use bootstrapping as a fallback.

if nargin<3,
    fcn = @mean;
end

nFiles = max(idx);

% Get the standard deviation across files
if nFiles>1,
    values = zeros( nFiles, 1 );

    for i=1:nFiles,
        values(i) = fcn( stat(idx==i) );
    end

    stdev = std(values);

% If only one file is provided, use bootstrap samples as a fallback.
% This is not comparable to standard deviation, but oh well.
else
    stdev = std( bootstrp(100,fcn,stat) );
end

end %function stdbyfile



function histdata = survival(times,dwellaxis)
% Make and display a survival plot for exponentially distributed
% dwell-times. dwellaxis are the bins to use for the histogram.

if numel(times)<1,
    histdata = [];
    return;
end

histdata = histc( times, dwellaxis );
histdata = sum(histdata) - cumsum(histdata);  %survival plot
histdata = histdata/histdata(1);  %normalize

plot( dwellaxis, histdata );
    
end



function output = expfit(times,dwellaxis)
% Y = a*exp(b*x)

histdata = histc( times, dwellaxis );
histdata = sum(histdata) - cumsum(histdata);  %survival plot
histdata = histdata/histdata(1);  %normalize
    
f = fit( dwellaxis',histdata, 'exp1', 'StartPoint',[1 -1/mean(times)] );
output = -1/f.b;

end


