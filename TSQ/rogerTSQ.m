% function rogerTSQ()
% TODO: add errors bars for stats, use stretched exponential fit for dwell times
% to account for heterogeneity. The stretch factor could be plotted also as an
% error to show the spread in the behavior. Beta closer to 0 would give high
% errors, close to 1 would give no error. (1-Beta)*tau


% PROCEEDURE:

% 1) Request filenames from user.
fnames = getFiles('*.txt');
if isempty(fnames), return; end
nFiles = numel(fnames);

%
initialModel = qub_loadModel;


%% % 2) Get the kinetic model used for analyzing the data...


% prep output variables
intensity    = zeros( nFiles,1 );
intensityStd = zeros( nFiles,1 );
SNRs         = zeros( nFiles,1 );
SNRsStd      = zeros( nFiles,1 );
tON       = zeros( nFiles,1 ); % dwell time of the first ON state dwell.
names     = cell(  nFiles,1 );
dwtFilename = cell(  nFiles,1 );

for i=1:nFiles,
    
    % 2) Load data from each file and calculate stats.
    % (Should data be filtered first??).
    % Donor intensity is scaled, so undo the scaling.
    constants = cascadeConstants;
    
    [d,a,f,~,time] = loadTraces( fnames{i} );
    data = d/constants.gamma+a;
    [nTraces,traceLen] = size(d);
    [~,names{i}] = fileparts( fnames{i} );
    
    stats = traceStat( d,a,f );
    
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
    
    % 3) Idealize the intensity data to a two-state model using SKM.
    sampling = ( time(2)-time(1) )*1000; %exposure time in seconds.
    
    skmParams.seperately = 1;
    skmParams.quiet = 1;
        
    [dwt,newModel,LL,offsets] = skm( data, sampling, initialModel, skmParams );
    dwtFilename{i} = strrep( fnames{i}, '.txt','.qub.dwt' );
    
%     mu    = newModel.mu(newModel.class);
%     sigma = newModel.sigma(newModel.class);
    mu    = initialModel.mu';
    sigma = initialModel.sigma';
    FRETmodel = [mu sigma]
    saveDWT( dwtFilename{i}, dwt, offsets, FRETmodel, sampling/1000 );
end

%%

% Calculate dwell times by fitting to exponential decays..
params.useCorrectedDwelltimes = 0; %don't merge blinks)
tON = lifetime_exp( dwtFilename, params );
tON = tON(:,2);

% Save results to a file that can be plotted as bar graphs in Origin.
fid = fopen('tsqStats.txt','w');
fwrite(fid, sprintf('Name\tIntensity (photons)\tIntensity stdev\tSNR\tSNR stdev\ttON (sec)\r\n') );

for i=1:nFiles,
    fwrite(   fid, sprintf('%s\t%.0f\t%.0f\t%.1f\t%.1f\t%.2f\r\n', ...
              names{i}, intensity(i),intensityStd(i), SNRs(i),SNRsStd(i), tON(i) )   );
end

fclose(fid);



% end %FUNCTION rogerTSQ




