function rtd_pipeline()

%% ANALYSIS PIPELINE FOR POST-SYNCHRONIZED tRNA SELECTION EXPERIMENTS
% To use this code, first go to the directory 

clear all;

constants = cascadeConstants;
sampling = 10; %ms


% Get the location of the data to process from the user.
% All data to be processed together should be in the same directory.
dataDir = uigetdir(pwd,'Select data directory');
cd(dataDir);

% Get the locations of kinetic models used for idealization.
modelDir = uigetdir(constants.modelLocation,'Select directory where tRNA selection models are located:');


% Load model and set re-estimation constraints: based 2-state model.
model1 = qub_loadModel( [modelDir filesep '090520_FretHist_2State_model.qmf'] );
model1.fixMu    = ones( model1.nStates,1 );
model1.fixSigma = ones( model1.nStates,1 );
fretModel1 = [model1.mu' model1.sigma'];

% Load model and set re-estimation constraints: full 6-state kinetic model
model2 = qub_loadModel( [modelDir filesep '090520_OH_initial_model_k20linear.qmf'] );
model2.fixMu    = ones( model2.nStates,1 );
model2.fixSigma = ones( model2.nStates,1 );
fretModel2 = [model2.mu' model2.sigma'];

% Parameters for SKM
skmParams.maxItr = 10;
skmParams.convLL = 0.01;
skmParams.zeroEnd = 1;
skmParams.quiet = 1;
skmParams.fixRates = 1;


%% 3. Run gettraces to extract traces from movies
disp('Running gettraces...');

gettracesParams.nPixelsToSum = 5; %old way
gettracesParams.don_thresh = 1800;
gettracesParams.overlap_thresh = 1.8;
gettracesParams.skipExisting = 1;
gettraces(dataDir,gettracesParams);


%% 4. Run autotrace to filter acquired data
disp('Running autotrace...');

% Find the location of all traces files to analyze
d = dir( [dataDir filesep '*.traces'] );
tracesFiles = {d.name};
nFiles = numel(tracesFiles);

% Define selection criteria.
selectionCriteria.min_maxFRET = 0.14;   %at least one frame with FRET>0.14
selectionCriteria.min_lifetime = 20;    %donor lifetime > 20 frames
selectionCriteria.min_snr = 6;          %SNR1 at least 6:1
selectionCriteria.max_ncross = 3;       %less than 4 donor blinking events
selectionCriteria.max_firstFRET = 0.14; %ignore accommodated traces.

% Select traces according to selection criteria
d_out = []; a_out = []; f_out = []; ids_out = {};

for i=1:nFiles,

    % Load traces file
    [d,a,f,ids,time] = loadTraces( tracesFiles{i} );
    
    % Select traces according to selection criteria.
    stats = traceStat( d,a,f );
    indexes = pickTraces( stats, selectionCriteria );
    
    d_out   = [d_out ; d(indexes,:)];
    a_out   = [a_out ; a(indexes,:)];
    f_out   = [f_out ; f(indexes,:)];
    ids_out = [ids_out ids(indexes)];

end

% Save selected traces
saveTraces( 'selected_traces.txt', 'txt', d_out,a_out,f_out,ids_out,time );



%% 5. Seperate events

disp('Seperating events...');

autosort( 'selected_traces.txt' );
drawnow;





%% ---- 6. Make a preliminary 1D post-synchronized histogram...

% Load FRET data (seperated events)
[donor,acceptor,fret,ids,time] = loadTraces('ips.txt');
[nTraces,traceLen] = size(fret);

% Idealize FRET data using SKM
[dwt] = skm( fret, sampling, model1, skmParams );

% Remove traces with no transitions
nDwells = cellfun( @numel, dwt )./2;
selected = nDwells>1;

dwt     = dwt( selected );
nTraces = sum(selected);
offsets = traceLen*((1:nTraces)-1);

donor    = donor( selected, :);
acceptor = acceptor( selected, :);
fret     = fret( selected, :);
ids      = ids( selected );

% Save the filtered data...
saveTraces( 'ips.flt.txt', 'txt', donor,acceptor,fret,ids,time );
forQuB2( {'ips.flt.txt'} );
saveDWT( 'ips1DHst.flt.qub.dwt', dwt, offsets, fretModel1, sampling );
%save selection list??

% Plot the results
frethist_norm_10ms( 'ips.flt.txt' );
title('Basic Post-synchronized Histogram');


% TODO: Create FRET histograms for Origin



%% Generate a list file
[nTraces,traceLen] = size(fret);
x = repmat( [1 traceLen-1],[nTraces,1] );
% y = repmat( offsets', [1,2] );
y = repmat( [0:traceLen:(nTraces*traceLen-1)]', [1,2] );
list = x+y-1;
dlmwrite( 'ips1DHst.flt.lst.txt',list, ...
          'delimiter','-', 'precision','%-10.f' );

% Make 1D histograms from dark and non-zero FRET states
statehist_2stateModel_v3( 'ips1DHst.flt.qub.dwt', ...
           'ips.flt.qub.txt', 'ips1DHst.flt.lst.txt', ...
           10, 420);


%% ---- 7. Identify outliers in QUB

% Idealize the filtered data using SKM
clear dwt; clear offsets;
[dwt] = skm( fret, sampling, model2, skmParams );

% Remove outliers
nDwells = cellfun( @numel, dwt )./2;
selected = nDwells>1 & nDwells<( mean(nDwells)+2*std(nDwells) );

% Create filtered dataset
dwt     = dwt( selected );
nTraces = sum(selected);
offsets = traceLen*((1:nTraces)-1);

donor    = donor( selected, :);
acceptor = acceptor( selected, :);
fret     = fret( selected, :);
ids      = ids( selected );

% Save the filtered data...
saveTraces( 'ips.flt.flt.txt', 'txt', donor,acceptor,fret,ids,time );
forQuB2( {'ips.flt.flt.txt'} );
saveDWT( 'ips.flt.flt.qub.dwt', dwt, offsets, fretModel2, sampling );




%% ---- 9. Cut traces after peptide bond formation time

% Seperate data into three groups, according to whether a peptide bond was
% formed or not: allMol12, Pep120, and noPep120.
cuttraces_2_10ms( 'ips.flt.flt.txt', 'ips.flt.flt.qub.dwt' );


% Re-idealize each of the datasets. This is neccessary because the data
% were modified by cuttraces. This also enables the model used for
% idealization to better fit the unique features of each dataset.
[d,a,fret] = loadTraces('allMol_ac120.txt');
[dwt,m,l,offsets] = skm( fret, sampling, model2, skmParams );
saveDWT( 'allMol_ac120.qub.dwt', dwt, offsets, fretModel2, sampling );
copyfile( 'allMol_ac120.qub.dwt', 'allMol_ac120-noise.qub.dwt' );

[d,a,fret] = loadTraces('PEP120.txt');
[dwt,m,l,offsets] = skm( fret, sampling, model2, skmParams );
saveDWT( 'PEP120.qub.dwt', dwt, offsets, fretModel2, sampling );

[d,a,fret] = loadTraces('noPep120.txt');
[dwt,m,l,offsets] = skm( fret, sampling, model2, skmParams );
saveDWT( 'noPep120.qub.dwt', dwt, offsets, fretModel2, sampling );


%% ---- 10. Make plots for viewing the results of idealization

% Make FRET contour plots for each dataset.
frethist_norm_10ms( 'allMol_ac120-noise.txt' );
title('Post-synchronized AC120 Histogram');

frethist_norm_10ms( 'PEP120.txt' );
title('Post-synchronized PEP120 Histogram');

frethist_norm_10ms( 'noPep120.txt' );
title('Post-synchronized NoPep120 Histogram');

%% Make transition-density plots for each dataset.
tdpOptions = {'normalize','total transitions','tdp_max',1.0 };
figure;
set(gcf,'Position',[166 391 1037 404]);

subplot(1,3,1);
tdpData = tdplot_rtdata_10ms( 'allMol_ac120.qub.dwt','allMol_ac120.txt' );
tplot( tdpData, tdpOptions{:} );
xlim([-0.15 0.8]); xlabel('Initial FRET');
ylim([-0.15 0.8]); ylabel('Final FRET');
title('All Molecules');

subplot(1,3,2);
tdpData = tdplot_rtdata_10ms('PEP120.qub.dwt','PEP120.txt');
tplot( tdpData, tdpOptions{:} );
xlim([-0.15 0.8]); xlabel('Initial FRET');
ylim([-0.15 0.8]); ylabel('Final FRET');
title('PEP120');

subplot(1,3,3);
tdpData = tdplot_rtdata_10ms('noPep120.qub.dwt','noPep120.txt');
tplot( tdpData, tdpOptions{:} );
xlim([-0.15 0.8]); xlabel('Initial FRET');
ylim([-0.15 0.8]); ylabel('Final FRET');
title('noPep120');



%%
end


