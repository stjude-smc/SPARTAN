function tsqComparison(files)



nBins = 40;  

% Get list of files to compare from user
if ~exist('files','var'),
    files = getFiles
end
nFiles = numel(files);


%     sampling = 40;
f = inputdlg('What is the sampling interval (in ms) for this data?');
sampling = str2double(f)


% 
% cumsum histograms of lifetime must have same axes for each file
% in terms of frames.  for now, just set aside space for number of bins.
names = cell(nFiles,1);
meanIntensity = zeros(nFiles,1);
donorLifetime = ones(nFiles,nBins+1);
fretLifetime  = ones(nFiles,nBins+1);

binCenters = [];

for i=1:nFiles,
   
    % Load dataset
    [d,a,f] = loadTraces(files{i});
    
    % Calculate trace properties, including intensity and bleaching time
    stats = traceStat(d,a,f);
    
    % 
    meanIntensity(i) = mean( [stats.t] );
    
    % donor lifetime
    LTdonor = [stats.lifetime] * (sampling/1000);
    
    if isempty(binCenters),
        [dlife,binCenters] = hist( LTdonor, nBins );
    else
        dlife = hist( LTdonor, binCenters );
    end
    
    dlife = cumsum(dlife)/sum(dlife);
    donorLifetime(i,2:end) = 1-dlife;
    
    % acceptor lifetime
    LTacc = [stats.acclife] * (sampling/1000);
    alife = hist( LTacc, binCenters );
    
    alife = cumsum(alife)/sum(alife);
    fretLifetime(i,2:end) = 1-alife;
    
end

bins = [0 binCenters];

% Save intensity information
save( 'meanIntensity.txt', 'meanIntensity', '-ASCII' );

% Save donor bleaching information
fl = [bins' donorLifetime'];
save( 'donorlife.txt', 'fl', '-ASCII' );

% Save acceptor bleaching information
fl = [bins' fretLifetime'];
save( 'acclife.txt', 'fl', '-ASCII' );




