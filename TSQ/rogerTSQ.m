function rogerTSQ(fnames)
% TODO: add errors bars for stats, use stretched exponential fit for dwell times
% to account for heterogeneity. The stretch factor could be plotted also as an
% error to show the spread in the behavior. Beta closer to 0 would give high
% errors, close to 1 would give no error. (1-Beta)*tau


% PROCEEDURE:


%-------------------------------------------------------------
% 1) Request filenames from user.
if nargin<1,
    fnames = getFiles;
    if isempty(fnames), return; end
end
nFiles = numel(fnames);

%
initialModel = qub_loadModel;
initialModel.fixSigma = [1 0];

constants = cascadeConstants;
constants.gamma=1;


% prep output variables
intensity    = zeros( nFiles,1 );
intensityStd = zeros( nFiles,1 );
SNRs         = zeros( nFiles,1 );
SNRsStd      = zeros( nFiles,1 );
%pt   = zeros(nFiles,1);
Ton  = zeros(nFiles,1);
% Toff = zeros(nFiles,1);
totalTon     = zeros( nFiles,1 );
%totalToff    = zeros( nFiles,1 );

names     = cell(  nFiles,1 );
dwtFilename = cell(  nFiles,1 );

fh = figure;

for i=1:nFiles,
    
    %-------------------------------------------------------------
    % 2) Load data from each file and calculate stats. 
    [d,a,f,ids,time] = loadTraces( fnames{i} );
    [~,names{i},ext] = fileparts( fnames{i} );
    
    
    %-------------------------------------------------------------
    % 3) Pick traces according to defined criteria if not already filtered.
    stats = traceStat( d,a,f, constants );
    
    if strcmp(ext,'.traces'),
        % Fit the background distribution to find a good cutoff.
        % This removes traces with low-intensity resurrection that is not
        % well-handled by the HMM analysis.
        bins = 0:20:2000;
        bg = [stats.bg];
        bg = bg(bg>0&bg<2000);
        N = hist( bg, bins );
        fitResult = fit( bins', N', 'gauss1', 'StartPoint',[25,median(bg),1.75*median(bg)] );
        cutoff = fitResult.b1 + 4*fitResult.c1; %4 standard deviations from mean.
        
        cla; hold on; bar(bins,N); plot(fitResult);
        plot( [cutoff cutoff], [0 max(N)], 'g-' );
        
        % Define selection criteria
        criteria.min_snr = 10;
        criteria.overlap = 1;
        criteria.max_bg = cutoff

        % Select traces according to criteria defined above.
        indexes = pickTraces( stats, criteria );
        d = d(indexes,:);
        a = a(indexes,:);
        f = f(indexes,:);
        ids = ids(indexes);
        stats = stats(indexes);
    end
    
    [nTraces,traceLen] = size(d);
    data = d+a;
    
    
    % Save stats by fitting the distributions.
    t = [stats.t];
%     bins = 0:500:30000;
%     [histdata] = hist( t(t>0), bins );
%     histdata = histdata / sum(histdata);
%     f = fit( bins',histdata', 'gauss1' );
%     rawIntensity = f.b1;
    rawIntensity = median(t);
    intensity(i) = rawIntensity*3.1/100; %adjusted for 100x gain and ADU/photon conversion.
    intensityStd(i) = std(t)*3.1/100;
    
%     if abs(rawIntensity-median(t))/rawIntensity > 0.15,
%         warning('Intensity fit is uncertain. May be multiple peaks!');
%     end
    
    snr = [stats.snr_s];
%     bins = 0:1:40;
%     [histdata] = hist( snr(snr>0), bins );
%     histdata = histdata / sum(histdata);
%     f = fit( bins',histdata', 'gauss1' );
%     SNRs(i) = f.b1; %mean
    SNRs(i) = median(snr);
    SNRsStd(i) = std(snr);
    
%     if abs(SNRs(i)-median(snr))/SNRs(i) > 0.15,
%         warning('SNR fit is uncertain. May be multiple peaks!');
%     end
    
    % Need to scale the data so it fits
    data = data/rawIntensity;
    
    
    %-------------------------------------------------------------
    % 4) Idealize the intensity data to a two-state model using SKM.
    sampling = ( time(2)-time(1) ); %exposure time in seconds.
    
    skmParams.seperately = 1;
    skmParams.quiet = 1;
        
    [dwt,newModel,LL,offsets] = skm( data, sampling, initialModel, skmParams );
        
    
    %-------------------------------------------------------------
    % 5) Remove traces that are poorly idealized or are significant outliers.
    % A) Remove any trace with >2 stdev transitions/total lifetime.
    % WARNING: assumes single-channel data!
    donorlife = [stats.donorlife];
    nTransitions = cellfun('size',dwt,1)-1;
    donorlife = reshape( donorlife, numel(donorlife),1 );
    nTransitions = reshape( nTransitions, numel(nTransitions),1 );
    
    transRate = nTransitions;%./donorlife;
    meanTrans = mean(transRate);
    stdTrans  = std(transRate);
    
    if stdTrans>0,
        selected = transRate < (meanTrans + 2*stdTrans);
        nRejected = nTraces-sum(selected);
    else
        selected = logical( ones(size(transRate)) );
        nRejected = 0;
    end
    
    disp( sprintf('%d) Removed %d traces (%.1f%%)', i, nRejected, 100*nRejected/nTraces) );
  
    
    %-------------------------------------------------------------
    % 6) Save idealization.
%     mu    = newModel.mu(newModel.class);
%     sigma = newModel.sigma(newModel.class);
    mu    = initialModel.mu';
    sigma = initialModel.sigma';
    FRETmodel = [mu sigma];
    dwtFilename{i} = [ removeExt(fnames{i}) '.qub.dwt' ];
    offsets = traceLen*((1:sum(selected))-1);
    saveTraces( [removeExt(fnames{i}) '_auto.txt'], 'txt', ...
                d(selected,:),a(selected,:),f(selected,:),ids(selected),time );
    saveDWT( dwtFilename{i}, dwt(selected), offsets, FRETmodel, sampling );
    
    offsets = traceLen*((1:sum(~selected))-1);
    saveTraces( [removeExt(fnames{i}) '_rejected.txt'], 'txt', ...
                d(~selected,:),a(~selected,:),f(~selected,:),ids(~selected),time );
    saveDWT( strrep(dwtFilename{i},'.dwt','_rejected.dwt'), dwt(~selected), ...
             offsets, FRETmodel, sampling );
         
         
    %-------------------------------------------------------------   
    % 6) Calculate total lifetimes
    
    % Figure out which class is the ON state
    means = FRETmodel(:,1);
    onState  = find( means==max(means) );
    %offState = find( means==min(means) );
    
    % Collect list of dwell-times in each class
    %nDwells = sum( cellfun(@numel,dwt ) )/2;
    onTimes = [];
    %offTimes = [];
    totalOnTimes  = zeros(nTraces,1);
    %totalOffTimes = zeros(nTraces,1);
    
    for j=1:nTraces,
        classes = dwt{j}(:,1);
        times   = dwt{j}(:,2).*sampling/1000;
        
        onTimes  = [onTimes ; times(classes==onState) ];
        %offTimes = [offTimes ; times(classes==offState) ];
        
        totalOnTimes(j)  = sum( times(classes==onState ) );
        %totalOffTimes(j) = sum( times(classes==offState) );
    end
    
    % Calculate percent time on
    %total_on = sum(onTimes);    
    %pt(i) = total_on / ( total_on + sum(offTimes) );
    
    
    % Calculate average total time on/off
    Ton(i)  = mean( onTimes  );
    %Toff(i) = mean( offTimes );    
    
    % Calculate average total time on/off
    totalOnTimes  = totalOnTimes(  totalOnTimes~=0  );
    %totalOffTimes = totalOffTimes( totalOffTimes~=0 );
    
    totalTon(i)  = mean( totalOnTimes  );
    %totalToff(i) = mean( totalOffTimes );
end

%%

% Calculate dwell times by fitting to exponential decays..
params.useCorrectedDwelltimes = 0; %don't merge blinks)
params.fitSingle=true;
Ton = lifetime_exp( dwtFilename, params );
Ton = Ton(:,2);

% Save results to a file that can be plotted as bar graphs in Origin.
f=0;k=0;
while all(f==0) && k<2,
    [f,p] = uiputfile('tsqStats.txt','Choose a filename to save the results.');
    fname = [p f];
    k = k+1;
end

% if fname,
    fid = fopen(fname,'w');

    fwrite(fid, sprintf('Name\tIntensity (photons)\tIntensity stdev\tSNR\tSNR stdev\ttON (sec)\tTotal tON (sec)\r\n') );

    for i=1:nFiles,
        fwrite(   fid, sprintf('%s\t%.0f\t%.0f\t%.1f\t%.1f\t%.2f\t%.2f\r\n', ...
                  names{i}, intensity(i),intensityStd(i), SNRs(i),SNRsStd(i), Ton(i), totalTon(i) )   );
    end
    
    fclose(fid);
% end

close(fh);


end
% end %FUNCTION rogerTSQ


function filename = removeExt( filename )

[p,n] = fileparts(filename);
filename = [p filesep n];

end

