function rtd_pipeline()

%% ANALYSIS PIPELINE FOR POST-SYNCHRONIZED tRNA SELECTION EXPERIMENTS
% To use this code, first go to the directory 


% Get the location of the data to process from the user.
% All data to be processed together should be in the same directory.
dataDir = uigetdir(pwd,'Select data directory');
cd(dataDir);

% modelDir = uigetdir(pwd,'Select directory where tRNA selection models are located:');
modelDir = '/home/dsterry/Documents/1_PROJECTS/RTD automation/090520 SKM test/models/';



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



%% ---- 6. Make a preliminary 1D post-synchronized histogram (OPTIONAL)

sampling = 10; %ms

% Load model and set re-estimation constraints
model1 = qub_loadModel( [modelDir filesep '090520_FretHist_2State_model.qmf'] );
model1.fixMu    = ones( model1.nStates,1 );
model1.fixSigma = ones( model1.nStates,1 );
fretModel1 = [model1.mu' model1.sigma']';

% Load FRET data (seperated events)
[donor,acceptor,fret,ids,time] = loadTraces('ips.txt');
assert( time(2)-time(1)==1 || time(2)-time(1)==10, 'Expecting 10ms data.' );

% Idealize FRET data using SKM
skmParams.maxItr = 100;
skmParams.convLL = 0.001;
skmParams.zeroEnd = 1; %this is critical for cuttraces().
[dwt,newModel,LL,offsets] = skm( fret, sampling, model1, skmParams );

% Remove traces with no transitions
nDwells = cellfun( @numel, dwt )./2;
selected = nDwells>1;

dwt      = dwt( selected );
offsets  = offsets( selected );

donor    = donor( selected, :);
acceptor = acceptor( selected, :);
fret     = fret( selected, :);
ids      = ids( selected );

% Save the filtered data...
saveTraces( 'ips.flt.txt', 'txt', donor,acceptor,fret,ids,time );
forQuB2( {'ips.flt.txt'} );
saveDWT( 'ips1DHst.flt.qub.dwt', dwt, offsets, fretModel1, sampling );

% Generate a list file for statehist.
[nTraces,traceLen] = size(fret);
x = repmat( [1 traceLen-1],[nTraces,1] );  %start and end points
y = repmat( [0:traceLen:(nTraces*traceLen-1)]', [1,2] );  %offsets
list = x+y-1;
dlmwrite( 'ips1DHst.flt.lst.txt',list, ...
          'delimiter','-', 'precision','%-10.f' );

% Make 1D histograms for dark and non-zero FRET states.
% This can then be used for fitting in Origin to define states.
statehist_2stateModel_v3( 'ips1DHst.flt.qub.dwt', ...
           'ips.flt.qub.txt', 'ips1DHst.flt.lst.txt', 10, 420);
       
% Make population contour plots (time-FRET histograms)
frethist_norm_10ms( 'ips.flt.txt' );
title('Basic Post-synchronized Histogram');
drawnow;
       



%% ---- 7. Identify outliers in QUB

% Load model and set re-estimation constraints
clear fretModel;
model2 = qub_loadModel( [modelDir filesep '090520_OH_initial_model_k20linear.qmf'] );
model2.fixMu    = ones( model2.nStates,1 );
model2.fixSigma = ones( model2.nStates,1 );

fretModel2 = [model2.mu' model2.sigma']';

% Idealize the filtered data using SKM
clear dwt; clear offsets;
[dwt,newModel,LL,offsets] = skm( fret, sampling, model2, skmParams );

% Remove outliers
nDwells = cellfun( @numel, dwt )./2;
selected = nDwells>1 & nDwells<( mean(nDwells)+2*std(nDwells) );

dwt      = dwt( selected );
offsets  = offsets( selected );

donor    = donor( selected, :);
acceptor = acceptor( selected, :);
fret     = fret( selected, :);
ids      = ids( selected );

% Save the filtered data...
saveTraces( 'ips.flt.flt.txt', 'txt', donor,acceptor,fret,ids,time );
forQuB2( {'ips.flt.flt.txt'} );
saveDWT( 'ips.flt.flt.qub.dwt', dwt, offsets, fretModel2, sampling );




%% ---- 9. Cut traces after peptide bond formation time
tic
cuttraces_2_10ms( 'ips.flt.flt.txt', 'ips.flt.flt.qub.dwt' );
forQuB2( {'allMol_ac120.txt','allMol_ac120-noise.txt','PEP120.txt','noPep120.txt'} );
disp(toc);

frethist_norm_10ms( 'allMol_ac120-noise.txt' );
title('Post-synchronized AC120 Histogram');
frethist_norm_10ms( 'PEP120.txt' );
title('Post-synchronized PEP120 Histogram');
frethist_norm_10ms( 'noPep120.txt' );
title('Post-synchronized NoPep120 Histogram');


%% ---- 10. Generate FRET histograms

frethist_norm_10ms( 'allMol_ac120-noise.txt' );
title('Post-synchronized AC120 Histogram');


% Load FRET data
[donor,acceptor,fret] = loadTraces('allMol_ac120-noise.txt');
[nTraces,traceLen] = size(fret);

% Idealize to a 2-state model so zero-FRET can be removed
[dwt,newModel,LL,offsets] = skm( fret, sampling, model1, skmParams );

% Save the results
saveDWT( 'allMol_ac120-noise.qub.dwt', dwt, offsets, fretModel1, sampling );

% Generate a list file
x = repmat( [1 traceLen-1],[nTraces,1] );
y = repmat( offsets', [1,2] );
list = x+y-1;
dlmwrite( 'allMol_ac120-noise.qub_sel.txt',list, ...
          'delimiter','-', 'precision','%-10.f' );

% Make 1D histograms from dark and non-zero FRET states
statehist_2stateModel_v3( 'allMol_ac120-noise.qub.dwt', ...
           'allMol_ac120-noise.qub.txt', 'allMol_ac120-noise.qub_sel.txt', ...
           10, 420);


%% ---- 11. Generate TD Plot for all molecules

% Idealize to full model
[dwt,newModel,LL,offsets] = skm( fret, sampling, model2, skmParams );
saveDWT( 'allMol120.qub.dwt', dwt, offsets, fretModel2, sampling );

% Generate TD plot for Origin
tdplot_rtdata_10ms('allMol120.qub.dwt','allMol_ac120-noise.qub.txt',...
                   'allMol_ac120-noise.qub_sel.txt');

% Display TD plot in MATLAB
figure;
tplot( tdplot('allMol120.qub.dwt','allMol_ac120-noise.txt') );
% tplot( load('allMol_ac120-noise_qub_tdp.txt') );
xlim([-0.15 0.8]); xlabel('Initial FRET');
ylim([-0.15 0.8]); ylabel('Final FRET');
set(gcf,'Position',[562 454 342 337]);



%% 12. transition density matrix (optional)

% transmat_uncor







