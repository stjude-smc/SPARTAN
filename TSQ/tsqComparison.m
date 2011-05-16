function [tBleachAcceptor,tBleachDonor] = tsqComparison(files,titles)



nBins = 40;  

% Get list of files to compare from user
if nargin<1,
    files = getFiles;
end
nFiles = numel(files);

% Make titles for files if none are provided
if nargin<2
    % Remove underscores (subscript)
    titles = strrep(files,'_',' ');

    % Strip off path, leaving just filename
    for i=1:nFiles,
        titles{i} = strip_path( titles{i} );
    end
end

if ~iscell(files),  files={files};  end
if ~iscell(titles), files={titles}; end


%% Load data and calculate photobleaching and intensity information.

% cumsum histograms of lifetime must have same axes for each file
% in terms of frames.  for now, just set aside space for number of bins.
names = cell(nFiles,1);
meanIntensity = zeros(nFiles,1);
donorLifetime = ones(nFiles,nBins+1);
fretLifetime  = ones(nFiles,nBins+1);

binCenters = [];

for i=1:nFiles,
   
    % Load dataset
    [d,a,f,ids,time] = loadTraces(files{i});
    
    if ~exist('sampling','var')
        if time(1)==1,
            framerate = inputdlg('What is the sampling interval (in ms) for this data?');
            sampling = str2double(framerate)
        else
            sampling = time(2)-time(1)
        end
    end

    
    % Calculate trace properties, including intensity and bleaching time
    stats = traceStat(d,a,f);
    
    % 
    meanIntensity(i) = mean( [stats.t] );
    
    % donor lifetime
    LTdonor = [stats.donorlife] * (sampling/1000);
    
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


%% Save data for import into Origin.

bins = [0 binCenters]';
donorLifetime = donorLifetime';
fretLifetime  = fretLifetime';

% Save intensity information
save( 'meanIntensity.txt', 'meanIntensity', '-ASCII' );

% Save donor bleaching information
fl = [bins donorLifetime];
save( 'donorlife.txt', 'fl', '-ASCII' );

% Save acceptor bleaching information
fl = [bins fretLifetime];
save( 'acclife.txt', 'fl', '-ASCII' );


%% Display results to user for immediate interpretation.

% Show Donor photobleaching raw data
subplot( 1,2,1 );
plot( bins, donorLifetime,'LineWidth',2 );
title( 'Donor lifetime' );
ylabel('Fraction photobleached');
xlabel('Time (sec)');

% Show Acceptor photobleaching raw data
subplot( 1,2,2 );
plot( bins, fretLifetime,'LineWidth',2 );
title(' Acceptor Lifetime' );
ylabel('Fraction photobleached');
xlabel('Time (sec)');
legend( titles );

% Fit bleaching rates to exponential decays
tBleachDonor    = zeros(nFiles,1);
tBleachAcceptor = zeros(nFiles,1);

for i=1:nFiles
    expFit1 = fit( bins,donorLifetime(:,i), 'exp1' );
    expFit2 = fit( bins,fretLifetime(:,i), 'exp1' );
    
    tBleachDonor(i)    = -1/expFit1.b;
    tBleachAcceptor(i) = -1/expFit2.b;
end



end %function tsqComparison



%%
function output = strip_path( filename )

pos = find( filename==filesep );
output = filename(pos(end)+1:end);

end
