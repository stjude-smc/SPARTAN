function [results,labels] = blinkingKinetics( filenames, titles )
% Provides preliminary photophysical kinetic parameter estimates
% of smFRET data.
%
% Preliminary analysis of TSQ-dye conjugate experiments.
% 10/17/08

%   Copyright 2008-2015 Cornell University All Rights Reserved.


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

% 
if any( isempty(strfind(filenames,'.dwt')) )
    warning('Not all files appear to be .DWTs');
end


%% Load the data and process the results
pt = zeros(nFiles,1);
Ton = zeros(nFiles,1);
Toff = zeros(nFiles,1);
totalTon = zeros(nFiles,1);
totalToff = zeros(nFiles,1);


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
    
    % Figure out which class is the ON state
    means = model(1:2:end);
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


