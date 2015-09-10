function [results,labels] = totalLifetime( filenames, titles )
% Provides preliminary photophysical kinetic parameter estimates
% of smFRET data.
%
% Preliminary analysis of TSQ-dye conjugate experiments.
% 10/17/08

%   Copyright 2008-2015 Cornell University All Rights Reserved.


if nargin<3,
    nBins = 6000;
end

% If no files specified, prompt user for them.
if nargin==0,
    filenames = cell(0,1);

    disp('Select traces files, hit cancel when finished');
    while 1,
        [datafile,datapath] = uigetfile({'*.dwt'},'Choose a DWT file:');
        if datafile==0, break; end  %user hit "cancel"

        filenames{end+1} = fullfile(datapath,datafile);
        disp( filenames{end} );
    end
end

nFiles = numel(filenames);

if nFiles == 0,
    disp('No files specified, exiting.');
    return;
end


% Generate plot titles 
if ~exist('titles','var'),
    % Strip off path, leaving just filename
    for i=1:nFiles,
        [p,f] = fileparts( filenames{i} );
        f = strrep(f,'_',' ');
        f = strrep(f,' auto','');
        titles{i} = f;
    end
end



%% Load the data and process the results
pt = zeros(nFiles,1);
Ton = zeros(nFiles,1);
Toff = zeros(nFiles,1);
totalTon = zeros(nFiles,1);
totalToff = zeros(nFiles,1);

Ton_hist = zeros(nBins,nFiles);
Toff_hist = zeros(nBins,nFiles);
totalTon_hist = zeros(nBins,nFiles);

for i=1:nFiles,

    % Load Dwell-time information
    [p,f] = fileparts( filenames{i} );
%     dwtfname = strrep(filenames{i},'.txt','.qub.dwt');
    dwtfname = filenames{i};
    if ~exist(dwtfname,'file'),
        error('Can''t find DWT file for %s',f);
    end
        
    [dwells,sampling,offsets,model] = loadDWT( dwtfname );
    nTraces = numel(dwells);
    
    if numel(model)~=4,
        error('Can only analyze 2-class data');
    end
    
    if ~exist('dwellaxis','var'),
        dwellaxis = (0:1:(nBins-1))*sampling;
    end
    
    % Figure out which class is the ON state
    means = model(1:2:end)
    onState  = find( means==max(means) );
    offState = find( means==min(means) );
    
    % Collect list of dwell-times in each class
    onTimes  = [];
    offTimes = [];
    totalOnTimes  = zeros(nTraces,1);
    totalOffTimes = zeros(nTraces,1);
    
    for j=1:nTraces,
        classes = dwells{j}(:,1);
        times   = dwells{j}(:,2).*sampling;
        
        onTimes  = [ onTimes  ; times(classes==onState ) ];
        offTimes = [ offTimes ; times(classes==offState) ];
        
        totalOnTimes(j)  = sum( times(classes==onState ) );
        totalOffTimes(j) = sum( times(classes==offState) );
    end
    
    % Calculate percent time on
    total_on = sum(onTimes);    
    pt(i) = total_on / ( total_on + sum(offTimes) );
    
    % Calculate average total time on/off
    Ton(i)  = mean( onTimes  );
    Toff(i) = mean( offTimes );    
    
    % Calculate average total time on/off
    totalOnTimes  = totalOnTimes(  totalOnTimes~=0  );
    totalOffTimes = totalOffTimes( totalOffTimes~=0 );
    
    totalTon(i)  = mean( totalOnTimes  );
    totalToff(i) = mean( totalOffTimes );
    
    % Create survival plots
    Ton_hist(:,i) = survivalHist( onTimes, dwellaxis );   
    Toff_hist(:,i) = survivalHist( offTimes, dwellaxis ); 
    totalTon_hist(:,i) = survivalHist( totalOnTimes, dwellaxis );    
end

results = [  [Ton Toff totalTon totalToff]./1000  100*pt  ];
labels = {'Ton','Toff','Total Ton','Total Toff','Time On (%)'};

% Ask user to save results to file
[f,p] = uiputfile('blinking_kinetics.txt', 'Save kinetics data as...');
outfile = [p f];

fid = fopen(outfile,'w');
fprintf( fid, '%s\t', labels{:} );
fprintf( fid, '\n' );

for i=1:nFiles
    fprintf( fid, '%s\t', titles{i} );
    fprintf( fid, '%f\t', results(i,:) );
    fprintf( fid, '\n' );
end

fclose(fid);


% Save survival plots
saveTimes('Ton.txt', titles, [dwellaxis' Ton_hist] );
saveTimes('Toff.txt', titles, [dwellaxis' Toff_hist] );
saveTimes('totalTon.txt', titles, [dwellaxis' totalTon_hist] );



end



function survival = survivalHist( times, dwellaxis )
% Creates a survival plot with <len> bins of <bin> (ms) width.

data = histc( times./1000, dwellaxis );
survival = sum(data) - cumsum(data);
survival = survival/survival(1);

end


function saveTimes( filename, titles, output )

nFiles = numel(titles);

fid = fopen(filename,'w');

fprintf(fid,'Time (sec)');
for i=1:nFiles
    fprintf( fid, '\t%s', titles{i} );
    fprintf( fid, '\n' );
end

for i=1:size(output,1),
    fprintf(fid,'%f\t',output(i,:));
    fprintf(fid,'\n');
end

fclose(fid);

end

