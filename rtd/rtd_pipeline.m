

%% 3. Run gettraces to extract traces from movies

gettraces;

%% 4. Run autotrace to filter acquired data

autotrace;

%% 5. Seperate events

autosort;




%% ---- 6. Make a preliminary 1D post-synchronized histogram...
clear all;

modelDir = '/home/dsterry/Documents/1_PROJECTS/RTD automation/090520 SKM test/models/';
% modelDir = uigetdir('','Select directory where tRNA selection models are located:');

sampling = 10; %ms

% Load model and set re-estimation constraints
model1 = qub_loadModel( [modelDir filesep '090520_FretHist_2State_model.qmf'] );
model1.fixMu    = ones( model1.nStates,1 );
model1.fixSigma = ones( model1.nStates,1 );

fretModel1 = [model1.mu' model1.sigma']';

% Load FRET data (seperated events)
[donor,acceptor,fret,ids,time] = loadTraces('ips.txt');
[nTraces,traceLen] = size(fret);

% Idealize FRET data using SKM
params.maxItr = 10;
params.convLL = 0.01;
params.zeroEnd = 1;
[dwt,newModel,LL,offsets] = skm( fret, sampling, model1, params );

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
%save selection list??

% Plot the results
frethist_norm_10ms( 'ips.flt.txt' );
title('Basic Post-synchronized Histogram');


% Create FRET histograms for Origin

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

% Load model and set re-estimation constraints
clear fretModel;
model2 = qub_loadModel( [modelDir filesep '090520_OH_initial_model_k20linear.qmf'] );
model2.fixMu    = ones( model2.nStates,1 );
model2.fixSigma = ones( model2.nStates,1 );

fretModel2 = [model2.mu' model2.sigma']';

% Idealize the filtered data using SKM
clear dwt; clear offsets;
[dwt,newModel,LL,offsets] = skm( fret, sampling, model2, params );

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
disp(toc);


frethist_norm_10ms( 'allMol_ac120-noise.txt' );
title('Post-synchronized AC120 Histogram');
frethist_norm_10ms( 'PEP120.txt' );
title('Post-synchronized PEP120 Histogram');
frethist_norm_10ms( 'noPep120.txt' );
title('Post-synchronized NoPep120 Histogram');


forQuB2( {'allMol_ac120.txt','allMol_ac120-noise.txt','PEP120.txt','noPep120.txt'} );


%% ---- 10. Generate FRET histograms

% Load FRET data
[donor,acceptor,fret,ids,time] = loadTraces('allMol_ac120-noise.txt');
[nTraces,traceLen] = size(fret);

% Idealize to a 2-state model so zero-FRET can be removed
[dwt,newModel,LL,offsets] = skm( fret, sampling, model1, params );

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


%---- 11. Generate TD Plot

% Idealize to full model
[dwt,newModel,LL,offsets] = skm( fret, sampling, model2, params );
saveDWT( 'allMol120.qub.dwt', dwt, offsets, fretModel2, sampling );

%%
% ---- Generate TD plot
% Generate TD plot for Origin
tdplot_rtdata_10ms('allMol120.qub.dwt','allMol_ac120-noise.qub.txt',...
                   'allMol_ac120-noise.qub_sel.txt');

%% Display TD plot in MATLAB
figure;
tplot( tdplot('allMol120.qub.dwt','allMol_ac120-noise.txt'));
% tplot( load('allMol_ac120-noise.qub_tdp.txt') );
xlim([-0.15 0.8]); xlabel('Initial FRET');
ylim([-0.15 0.8]); ylabel('Final FRET');
set(gcf,'Position',[562 454 342 337]);



%% 12. transition density matrix (optional)

transmat_uncor;













