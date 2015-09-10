function dwt = hammyToDWT( datapath, sampling )
% Produces a DWT file just barely correct enough
% to produce TD plots...

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Get location of hammy analysis results (path.dat files) here
if nargin<1,
    datapath=uigetdir;
    if datapath==0, return; end
end

if nargin<2,
    f = inputdlg('What is the sampling interval (in ms) for this data?');
    sampling = str2double(f)
end

% Get list of files
pathFiles = dir( [datapath filesep '*path.dat'] );
nFiles = numel(pathFiles);

% Load the idealizations...
dwt = cell(nFiles,1);
h = waitbar(0,'Loading HaMMy results files...');

for i=1:nFiles
    % Load viterbi path
    filename = [datapath filesep pathFiles(i).name];
    data = load( [datapath filesep pathFiles(i).name] );
    
    if size(data,2)~=5
        warning('hammyToDWT:load','Malformed or empty file: %s',filename);
    end
    
    % Convert idealization to state trajectory
    fretPath = RLEncode( data(:,5)' );
    fret  = fretPath(:,1);
    fret( fret<0.1 ) = 0;  %dark state
    s = unique(sort(fret));
    
    for j=1:numel(s),
        fret( fret==s(j) ) = j;
    end
    
    times = fretPath(:,2);
    dwt{i} = [fret times];
    
    totalTime = sum(times);
    
    waitbar(i/nFiles);
end
close(h);

fakeModel = [0 0.24 0.39 0.56 0.83 ; 0.1 0.1 0.1 0.1 0.1]';

offsets = (0:(nFiles-1))*totalTime;

[f,p] = uiputfile('*.qub.dwt');
if f==0,
    [f,p] = uiputfile('*.qub.dwt');
end
saveDWT( [p f], dwt, offsets, fakeModel, sampling );
