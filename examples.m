% SPARTAN SCRIPTING EXAMPLES
% ==========================
% 
% This file contains examples for how to script SPARTAN functions.
% Each cell in this notebook can be run independently.
% Use these as templates for creating your own analysis routines.
%
% To run each example, click the "Run" button in the "Editor" tab of the
% MATLAB toolbar or Control+Enter on your keyboard with the desired example
% cell selected.
% 



%% Load and view raw movie data using MovieViewer
clear;

% Dialog to select input movie file (empty if user hits "Cancel")
filename = getFile('*.tif','Select a movie');

% Load file as as a Movie object. No image data is loaded yet.
mov = Movie.load(filename);

% Display a viewer UI with scroll bars for intensity and frame number
mv = MovieViewer(filename);
mv.show();


% Average the first 10 images in the stack. "frames" is a 3-dimensional
% matrix with dimensions of row, column, frame number.
frames = mov.readFrames(1:10);
stktop = mean(frames,3);

% Display image and line profiles
figure;
range = quantile(stktop(:), [0.1 0.9]);
imshow(stktop,range, 'Colormap',colormap('hot'));

% Display X and Y line profiles
figure;
subplot(1,2,1);
plot(1:mov.nX, sum(stktop,1));
xlabel('X Position');

subplot(1,2,2);
plot(1:mov.nY, sum(stktop,2));
xlabel('Y Position');



%% Isolating spectral channels using ChannelExtractor
clear;

filename = getFile('*.tif','Select a movie');

% Specifies how to split each image. The numbers define channel indexes and
% the shape of the matrix defines how to split it.
% [1] = single channel
% [2 3; 1 4] = four colors
% [1 2; 0 3] = three colors (ignore empty bottom right area)
fieldArrangement = [1 2];

% "channels" is a struct array describing each spectral band in the order
% defined by the indexes in "fieldArrangement". 
channels = struct('name',{'Cy3','Cy5'}, 'wavelength',{532,640}, 'photonsPerCount',0.49);

% Create an interface for extracting frame data with isolated channels.
chext = ChannelExtractor(filename, fieldArrangement, channels);
chext.skipFrames = 1;  %optional parameters are set like this

% Display the split movie in a UI with scroll bars.
viewer = MovieViewer(chext);
viewer.show();


% Average first 10 frames from each channel (donor and acceptor).
% "frames" is a cell array (one per channel) with the image data.
frames = chext.read(1:10);
donor = mean( frames{1}, 3);
acceptor = mean( frames{2}, 3);

% Create and overlay intensity line profiles
figure;
plot(1:chext.nX, sum(donor,1), 'g-');
hold on;
plot(1:chext.nX, sum(acceptor,1), 'r-');
xlabel('X Position');
legend({'Donor','Acceptor'});



%% Extract traces from movie data (MovieParser)
% This replicates a simple version of the "gettraces" program.
clear;

% Get the default instrument profile, which includes all the analysis
% parameters for the process. See cascadeConstants.m for details.
params = cascadeConstants('gettracesDefaultParams');
params.nAvgFrames = 20;  %adjust any parameters as needed like this

% Load movie and launch a UI for viewing frame data
filename = getFile('*.tif','Select a movie file');
mp = MovieParser(filename, params);
viewer = MovieViewer( mp.chExtractor );
viewer.show();

% Select particles as local maxima above automatically chosen threshold
mp.getPeaks();

% Display some summary statistics on the selected particles.
nPicked = size(mp.total_peaks,1);
nRejected = size(mp.rejectedTotalPicks,1);

fprintf('\nSelected %d/%d traces (%.1f%%)\n', nPicked, nRejected, ...
                                    100*nPicked/(nPicked+nRejected));
fprintf('Average PSF size: %.1f px\n',mp.psfWidth);
fprintf('Integration efficiency: %.1f%%\n',mp.integrationEfficiency);
fprintf('Alignment quality %.2f, deviation %.2f px.\n\n', ...
        mp.alignStatus(2).quality, mp.alignStatus(2).abs_dev);

% Graphically show peak centers
viewer.highlightPeaks( mp.peaks, mp.rejectedPicks, mp.total_peaks, mp.rejectedTotalPicks );

% Integrate traces and save as a .rawtraces file.
[p,f] = fileparts(filename);
outname = fullfile(p,[f '.rawtraces']);
mp.integrateAndSave(outname);



%% Select traces in one file according to defined criteria (autotrace)
% This script is equivalent to using autotrace with a single input file.
clear;

% Dialog to select input .traces file
filename = getFile('*.traces','Select input files');

% Load .traces file, creating a Traces object called 'data'.
data = loadTraces(filename);

% Calculate standard descriptive statistics from trace data.
% "stats" is a structure array, one element per trace, with fields named
% for each of the statistics (for example: "snr").
stats = traceStat(data);

% Define a set of criteria for selecting traces.
% Each field is constructed by combining a comparative operator and the
% statistic name from the "stats" struct above.
criteria.min_snr = 8;       % SNR_background > 8.
criteria.min_acclife = 50;  % at least 50 frames of FRET.
criteria.max_corr = 0.5;    % D/A correlation coefficient < 0.5.

% Select traces that meet all of the above criteria.
% The data object is modified "in place" (no copy is made).
% "selection" lists the indexes of traces that fit the criteria.
selection = pickTraces(stats,criteria);
data.subset(selection);

% It is also possible to create a copy named data2 instead like this: 
%data2 = data.getSubset(selection);

% Select traces by spatial location.
% Here, 'selection' is a logical array.
x = [data.traceMetadata.donor_x];
y = [data.traceMetadata.donor_y];
nX = data.fileMetadata.nX;
nY = data.fileMetadata.nY;

selection = x>nX*1/4 & x<nX*3/4 & y>nY*1/4 & y<nY*3/4;
data.subset(selection);

% Save this new trace data with a new file extension
[p,f,e] = fileparts(filename);
outname = fullfile(p,[f '_auto' e]);
data.save(outname);



%% Batch processing for trace selection (autotrace)
clear;

% Prompt user for input. "files" is a cell array of file paths.
files = getFiles('*.rawtraces');

% Uncomment to use the default criteria in autotrace
% criteria = cascadeConstants('defaultAutotraceCriteria');

% Manually define selection criteria instead
criteria.eq_overlap  = 0;       % Remove overlapping molecules
criteria.min_snr     = 8;       % SNR over background
criteria.min_acclife = 15;      % FRET lifetime

for i=1:numel(files)
    % Load traces and select a subset according to above criteria.
    data = loadTraces(files{i});
    selected = pickTraces( traceStat(data), criteria );
    data.subset(selected);

    % Save selected traces
    [p,f] = fileparts(files{i});
    outname = fullfile(p,[f '_auto.traces']);
    data.save(outname);
end



%% Hidden Markov Modeling: idealization and saving to .dwt (batchKinetics)
clear;

% Specify the filenames manually
traceFile = 'sim.traces';
modelFile = 'test.model';

% ...or use getFile() to prompt the user instead.
% modelFile = getFile('*.model','Select a model file');
% traceFile = getFile('*.traces','Select a traces file');

% Load a model file
model = QubModel(modelFile);

for i=1:model.nClasses
    fprintf('State %d: FRET = %.2f +/- %.2f\n', i, model.mu(i), model.sigma(i))
end

% Load data to analyze
data = loadTraces(traceFile);

% Idealize traces using the segmental k-means algorithm (SKM).
% 'idl' is a matrix the same size as the fret matrix (traces in rows, time
% points in columns) with the state assignment at each point. A value of
% zero means there is no idealization at that point.
% 'params' is optional.
idl = skm(data.fret, data.sampling, model);

% Save idealization to file
[p,f] = fileparts(traceFile);
outname = fullfile(p,[f '.dwt']);
[dwt,offsets] = idlToDwt(idl);
saveDWT(outname, dwt, offsets, model, data.sampling);


% You can display a model like this. Any changes made by the user in the
% viewer window will also modify the 'model' object in real time.
% viewer = QubModelViewer(model);



%% Assemble various plots into one figure (like makeplots)
clear;
figure;

% Modify to the path of your data...
traceFile = 'sim.traces';
%dwtFile = 'sim.dwt';

% Automatically find corresponding .dwt file, if it exists.
dwtFile = findDwt(traceFile);
dwtFile = dwtFile{1};
hasDwt = ~isempty(dwtFile);

rows = hasDwt+1;
cols = 4;

% Mean FRET efficiency over time
ax(1) = subplot(rows,cols,1);
avgFretTime(ax(1), traceFile);

% FRET histogram (1D)
ax(2) = subplot(rows,cols,2);
frethistComparison(ax(2), traceFile);

% FRET-time contour plot
ax(3) = subplot(rows,cols,3);
cplot(ax(3), traceFile);

% FRET-time contour plot with longer time axis than default
ax(4) = subplot(rows,cols,4);
options = cascadeConstants('defaultMakeplotsOptions');  %default settings
options.contour_length = 300;
cplot(ax(4), traceFile, options);
title('Long time axis');


% Plots that require idealized data.
% If the first argument (ax) is removed, each of these functions will
% launch their own window instead.
if hasDwt
    % State occupancy over time.
    % Will automatically search for associated .dwt file
    ax(4) = subplot(rows,cols,5);
    occtime(ax(4), traceFile);

    % FRET histogram for each state, overlayed.
    ax(5) = subplot(rows,cols,6);
    statehist(ax(5), dwtFile, traceFile);

    % Transition density plot
    ax(6) = subplot(rows,cols,7);
    tdp = tdplot(dwtFile, traceFile);
    tplot(ax(6), tdp);

    % Dwell-time histograms
    ax(7) = subplot(rows,cols,8);
    dwellhist(ax(7), dwtFile);
end



%% Overlaying plots from a list of files

% Get a cell array of file paths from the user
files = getFiles('*.traces');

% Side-by-side contour plots and histograms
makeplots(files);

% FRET histograms (1D) overlay
frethistComparison(files);

% Average FRET value over time
avgFretTime(files);

% Mean FRET value
meanFret(files);

% Fluorescence correlation plot
fluorcorr(files);


% ----------------------------
% These require idealized data

% Average transition rate
transitionsPerSecond(files);

% Average occupancy (overall fraction of time) in each state.
percentTime(files);

% dwell time histograms
dwellhist(files);



%% Simulating single molecule data
clear;

% Load a model from file
% model = QubModel('test.model');

% ...or create a new model
nStates = 3;
model = QubModel(nStates);
model.mu    = [0 0.4 0.7];       %mean fret values
model.sigma = [0.06 0.06 0.06];  %standard deviation of fret
model.p0    = [0 0.5 0.5];       %initial state probabilities
model.rates = [0    5   5;       %rate constants (row=from, col=to)
               0.1  0   1;
               0.1  1   0];
model.save('test.model');  %save to file for future reference

% Simulation parameters
nTraces = 1000;
nFrames = 3000;
sampling = 40;  %ms

options = struct('snr',30, 'shotNoise',true, ...
                 'gamma',1, 'crosstalk',0, 'ade',0, ...
                 'totalIntensity',500, 'stdTotalIntensity',0, ...
                 'accLife',20, 'donLife',30, 'simExpBleach',true );

% Run simulation and save results to file.
data = simulate(nTraces, nFrames, sampling/1000, model, options);
data.save('sim.traces');


% Display simulated traces in a new window
figure;
viewer = TraceListViewer(uipanel);
viewer.loadTraceData('sim.traces');
% view.loadIdealization('sim.dwt');





