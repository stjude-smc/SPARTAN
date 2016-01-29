function [tBleachAcceptor,tBleachDonor] = tsqComparison(files,titles)
%

%   Copyright 2007-2015 Cornell University All Rights Reserved.


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
        [~,titles{i}] = fileparts( titles{i} );
    end
end

if ~iscell(files),  files={files};  end
if ~iscell(titles), files={titles}; end

if nFiles==0, return; end


%% Load data and calculate photobleaching and intensity information.

% cumsum histograms of lifetime must have same axes for each file
% in terms of frames.  for now, just set aside space for number of bins.
meanIntensity = zeros(nFiles,1);
donorLifetime = ones(nFiles,nBins+1);
fretLifetime  = ones(nFiles,nBins+1);

binCenters = [];

for i=1:nFiles,
   
    % Load dataset
    data = loadTraces(files{i});
    
    if ~exist('sampling','var')
        if data.time(1)==1,
            framerate = inputdlg('What is the sampling interval (in ms) for this data?');
            sampling = str2double(framerate);
        else
            sampling = data.time(2)-data.time(1);
        end
    end

    
    % Calculate trace properties, including intensity and bleaching time
    stats = traceStat( data );
    clear data;
    
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

hf = figure;

% Show Donor photobleaching raw data
ax = subplot( 1,2,1, 'Parent',hf );
plot( ax, bins, donorLifetime,'LineWidth',2 );
title( ax, 'Donor lifetime' );
ylabel(ax, 'Fraction photobleached');
xlabel(ax, 'Time (s)');

% Show Acceptor photobleaching raw data
ax = subplot( 1,2,2, 'Parent',hf );
plot( ax, bins, fretLifetime,'LineWidth',2 );
title(ax, ' Acceptor Lifetime' );
ylabel(ax, 'Fraction photobleached');
xlabel(ax, 'Time (s)');
legend(ax, titles );

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



