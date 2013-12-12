function rtd_pipeline()

%% 0. ANALYSIS PIPELINE FOR POST-SYNCHRONIZED tRNA SELECTION EXPERIMENTS

%% 1. SET PIPELINE PARAMETERS

% Get the location of the data to process from the user.
% All data to be processed together should be in the same directory.
clear all;
constants = cascadeConstants;
dataDir = uigetdir(pwd,'Select data directory');
cd(dataDir);
disp('01. Set Pipeline Parameters...');
% Load model and set re-estimation constraints: based 2-state model.
model1 = qub_loadModel( [constants.modelLocation filesep 'tRNA selection/090520_FretHist_2State_model.qmf'] );
model1.fixMu    = ones( model1.nStates,1 );
model1.fixSigma = ones( model1.nStates,1 );
fretModel1 = [model1.mu' model1.sigma'];

% Load model and set re-estimation constraints: full 6-state kinetic model
model2 = qub_loadModel( [constants.modelLocation filesep 'tRNA selection/090520_OH_initial_model_k20linear.qmf'] );
model2.fixMu    = ones( model2.nStates,1 );
model2.fixSigma = ones( model2.nStates,1 );
fretModel2 = [model2.mu' model2.sigma'];

% Step 2 Parameters for gettraces
gettracesParams.nPixelsToSum = 4; %old way
gettracesParams.don_thresh = 1800;
gettracesParams.overlap_thresh = 1.8;
gettracesParams.skipExisting = 1;
gettracesParams.recursive = 1;
gettraces(dataDir,gettracesParams);

% Step3 Define selection criteria.
selectionCriteria.min_maxFret   = 0.14;  %at least one frame with FRET>0.14
selectionCriteria.min_lifetime  = 20;    %donor lifetime > 20 frames
selectionCriteria.min_snr       = 6;     %SNR1 at least 6:1
selectionCriteria.max_ncross    = 3;     %less than 4 donor blinking events
selectionCriteria.max_firstFret = 0.14;  %no post-accommodation traces.
%selectionCriteria.min_maxFret  = 0.55;  %must accommodate.
%selectionCriteria.min_acclife  = 100;

% Step 5: Parameters for SKM
skmParams.maxItr = 10;
skmParams.convLL = 0.01;
skmParams.zeroEnd = 1;
skmParams.seperately = 1;
skmParams.quiet = 1;
skmParams.fixKinetics = 1;

%Step 6: Parameters for Statehist
IntervalStatehist=32; %Integrate FRET-histogram over n images

%Step 8: Parameters for cuttraces
FRETcr=0.206; %FRET value of the CR-state
FRETga=0.335; %FRET value of the GA-state
FRETac=0.559; %FRET value of the AC-state
time2peptide=120; %unit: ms

%Step 9: Parameters for Fret-histogram
backset=4; %unit: Images
timeWindow=42; %unit: Images

%% 2. Run gettraces.m to extract traces from movies
% All data contained in the selected folder *and subfolders* is analyzed.

disp('02. Running gettraces...');

% Get the sampling rate
d = rdir( [dataDir filesep '**' filesep '*.rawtraces'] );
tracesFiles = {d.name};
data = loadTraces( tracesFiles{1} );
sampling = data.time(2); %ms
disp(['sampling [ms]:', num2str(sampling)]);

%% 3. Select traces according to selection criteria
disp('03. Running autotrace...');
d = rdir( [dataDir filesep '**' filesep '*.rawtraces'] );

opt.outFilename = 'selected_traces.traces';
loadPickSaveTraces( {d.name}, selectionCriteria, opt );


%% 4. Run autosort.m and seperate events

disp('04. Running autosort...');
autosort( 'selected_traces.traces' );
drawnow;

%% 5. Remove traces with no transitions and make a preliminary 2D post-synchronized FRET histogram

disp('05. Running filter_1...');
% Load FRET data (seperated events)
data = loadTraces('ips.traces');
[nTraces,traceLen] = size(data.fret);

% Idealize FRET data using SKM
[dwt] = skm( data.fret, sampling, model1, skmParams );

% Remove traces with no transitions
nDwells = cellfun( @numel, dwt )./2;
selected = nDwells>1;

dwt     = dwt( selected );
nTraces = sum(selected);
offsets = traceLen*((1:nTraces)-1);

data.donor    = data.donor( selected, :);
data.acceptor = data.acceptor( selected, :);
data.fret     = data.fret( selected, :);
data.traceMetadata = data.traceMetadata( selected );

% Save the filtered data...
saveTraces( 'ips.flt.traces', 'traces', data );
forQuB2( {'ips.flt.traces'} );
saveDWT( 'ips1DHst.flt.qub.dwt', dwt, offsets, fretModel1, sampling );

% Plot the results
frethist_norm_v4('ips.flt.traces',1,0,sampling);


%% 6. Run statehist_2StateModel_v3.m and generate a 1D FRET histogram

disp('06. Running statehist...');
x = repmat( [1 traceLen-1],[nTraces,1] );
% y = repmat( offsets', [1,2] );
y = repmat( [0:traceLen:(nTraces*traceLen-1)]', [1,2] );
list = x+y-1;
dlmwrite( 'ips1DHst.flt.lst.txt',list, ...
          'delimiter','-', 'precision','%-10.f' );

%Parameters for statehist in ms
lwLimit=sampling;
upLimit=IntervalStatehist*sampling;
disp(['lwLimit [ms]:', num2str(lwLimit),'   upLimit [ms]:', num2str(upLimit)]);
% Make 1D histograms from dark and non-zero FRET states
statehist_2stateModel_v3( 'ips1DHst.flt.qub.dwt', ...
'ips.flt.qub.txt', 'ips1DHst.flt.lst.txt',lwLimit, upLimit,sampling,traceLen);


%% 7. Identify outliers in QUB

disp('07. Running filter_2...');
% Idealize the filtered data using SKM
clear dwt; clear offsets;
[dwt] = skm( data.fret, sampling, model2, skmParams );

% Remove outliers
nDwells = cellfun( @numel, dwt )./2;
selected = nDwells>1 & nDwells<( mean(nDwells)+2*std(nDwells) );

% Create filtered dataset
dwt     = dwt( selected );
nTraces = sum(selected);
offsets = traceLen*((1:nTraces)-1);

data.donor    = data.donor( selected, :);
data.acceptor = data.acceptor( selected, :);
data.fret     = data.fret( selected, :);
data.traceMetadata = data.traceMetadata( selected );

% Save the filtered data...
saveTraces( 'ips.flt.flt.traces', 'traces', data );
forQuB2( {'ips.flt.flt.traces'} );
saveDWT( 'ips.flt.flt.qub.dwt', dwt, offsets, fretModel2, sampling );


%% 8. Cut traces after peptide bond formation time

disp('08. Running cuttraces...');
% Seperate data into three groups, according to whether a peptide bond was
% formed or not: allMol12, Pep120, and noPep120.
cuttraces_v3('ips.flt.flt.traces', 'ips.flt.flt.qub.dwt',FRETac,time2peptide,sampling);

% Re-idealize each of the datasets. This is neccessary because the data
% were modified by cuttraces. This also enables the model used for
% idealization to better fit the unique features of each dataset.
data = loadTraces('allMol_ac120.traces');
fret = data.fret;
if (size(fret,1))>0
    [dwt,m,l,offsets] = skm( fret, sampling, model2, skmParams );
    saveDWT( 'allMol_ac120.qub.dwt', dwt, offsets, fretModel2, sampling );
    file_allMol_ac120=1;
else
    disp('file is emty');
    file_allMol_ac120=0;
end

data = loadTraces('PEP120.traces');
fret = data.fret;
if (size(fret,1))>0
    [dwt,m,l,offsets] = skm( fret, sampling, model2, skmParams );
    saveDWT( 'PEP120.qub.dwt', dwt, offsets, fretModel2, sampling );
    file_PEP120=1;
else
    disp('file is emty');
    file_PEP120=0;
end

data = loadTraces('noPEP120.traces');
fret = data.fret;
if (size(fret,1))>0
    [dwt,m,l,offsets] = skm( fret, sampling, model2, skmParams );
    saveDWT( 'noPEP120.qub.dwt', dwt, offsets, fretModel2, sampling );
    file_noPEP120=1;
else
    disp('file is emty');
    file_noPEP120=0;
end
disp(['time2Peptide [ms]:', num2str(time2peptide)]);
%% 9. Make plots for viewing the results of idealization

disp('09. Making FRET-Histograms...');
CONTOUR_LENGTH=50;
figure;
set(gcf,'Position',[186 411 1037 404]);

% Make FRET contour plots for each dataset.
if file_allMol_ac120>0
    subplot(1,3,1);
    [nAllMol,normfactor,frethst]=frethist_norm_v4('ips.flt.flt.traces',1,0,sampling);
    time_axis=frethst(1,2:end);
    fret_axis=frethst(2:end,1);
    con = 0:0.01:0.13;
    cmap=dlmread('frethist_colormap_peter.txt')/255;
    [C,hand]=contourf(time_axis(1:CONTOUR_LENGTH),fret_axis,frethst(2:end,2:CONTOUR_LENGTH+1),con);
    set(gca,'PlotBoxAspectRatio',[1.5 2 1]);
    set(hand,'LineColor','none');
    colormap(cmap);
    axis([-backset*(sampling/1000) timeWindow*(sampling/1000) -0.15 0.9]);
    xlabel('Time [s]');
    ylabel('FRET');
    zoom on;
    title('FRET Histogram (All Molecules)');
    str1='N =  ';
    str2=num2str(nAllMol);
    str=strcat(str1,str2);
    annotation('textbox',[0.265 0.815 0.074 0.05],'BackgroundColor',[1,1,1],...
    'string',str,'FontSize',8);
end


if file_PEP120>0
    subplot(1,3,2);
    [nPEP,~,frethst]=frethist_norm_v4('PEP120org.traces',2,0,sampling,normfactor);
    time_axis=frethst(1,2:end);
    fret_axis=frethst(2:end,1);
    con = 0:0.01:0.13;
    cmap=dlmread('frethist_colormap_peter.txt')/255;
    [C,hand]=contourf(time_axis(1:CONTOUR_LENGTH),fret_axis,...
    frethst(2:end,2:CONTOUR_LENGTH+1),con);
    set(gca,'PlotBoxAspectRatio',[1.5 2 1]);
    set(hand,'LineColor','none');
    colormap(cmap);
    axis([-backset*(sampling/1000) timeWindow*(sampling/1000) -0.15 0.9]);
    xlabel('Time [s]');
    ylabel('FRET');
    zoom on;
    title('FRET Histogram (PEP Events)');
    str1='N =  ';
    str2=num2str(nPEP);
    str=strcat(str1,str2);
    annotation('textbox',[0.545 0.815 0.074 0.05],'BackgroundColor',[1,1,1],...
    'string',str,'FontSize',8);
end

if file_PEP120>0
    subplot(1,3,3);
    [nNoPep,~,frethst]=frethist_norm_v4('noPEP120.traces',2,0,sampling,normfactor);
    time_axis=frethst(1,2:end);
    fret_axis=frethst(2:end,1);
    con = 0:0.01:0.13;
    cmap=dlmread('frethist_colormap_peter.txt')/255;
    [C,hand]=contourf(time_axis(1:CONTOUR_LENGTH),fret_axis,...
    frethst(2:end,2:CONTOUR_LENGTH+1),con);
    set(gca,'PlotBoxAspectRatio',[1.5 2 1]);
    set(hand,'LineColor','none');
    colormap(cmap);
    axis([-backset*(sampling/1000) timeWindow*(sampling/1000) -0.15 0.9]);
    xlabel('Time [s]');
    ylabel('FRET');
    zoom on;
    title('FRET Histogram (NoPEP Events)');
    str1='N =  ';
    str2=num2str(nNoPep);
    str=strcat(str1,str2);
    annotation('textbox',[0.825 0.815 0.074 0.05],'BackgroundColor',[1,1,1],...
    'string',str,'FontSize',8);
end

%% 10. Make transition-density plots for each dataset.
disp('10. Making TD-Plots...');
tdpOptions = {'normalize','total transitions','tdp_max',1.0 };
figure;
set(gcf,'Position',[166 391 1037 404]);

if file_allMol_ac120>0 
    subplot(1,3,1);
    [nTraces1, nTrans1,tdpData] = tdplot_rtdata_v3( 'allMol_ac120.qub.dwt','allMol_ac120.traces',1);
    tplot( tdpData, tdpOptions{:} );
    xlim([-0.15 1.0]); xlabel('Initial FRET');
    ylim([-0.15 1.0]); ylabel('Final FRET');
    title('All Molecules');
    str1='T =  ';
    str2=num2str(nTrans1);
    str3='T/N =  ';
    str4=num2str(nTrans1/nTraces1);
    strA=strcat(str1,str2);
    strB=strcat(str3,str4);
    annotation('textbox',[0.265 0.730 0.074 0.05],'BackgroundColor',[1,1,1],...
    'string',strA,'FontSize',8);
    annotation('textbox',[0.265 0.675 0.074 0.05],'BackgroundColor',[1,1,1],...
    'string',strB,'FontSize',8);
end

if file_PEP120>0
    subplot(1,3,2);
    [nTraces2, nTrans2,tdpData] = tdplot_rtdata_v3('PEP120.qub.dwt','PEP120.traces',2,nTrans1);
    tplot( tdpData, tdpOptions{:} );
    xlim([-0.15 1.0]); xlabel('Initial FRET');
    ylim([-0.15 1.0]); ylabel('Final FRET');
    title('PEP Events');
    str1='T =  ';
    str2=num2str(nTrans2);
    str3='T/N =  ';
    str4=num2str(nTrans2/nTraces2);
    strA=strcat(str1,str2);
    strB=strcat(str3,str4);
    annotation('textbox',[0.545 0.730 0.074 0.05],'BackgroundColor',[1,1,1],...
    'string',strA,'FontSize',8);
    annotation('textbox',[0.545 0.675 0.074 0.05],'BackgroundColor',[1,1,1],...
    'string',strB,'FontSize',8);
end

if file_noPEP120>0
    subplot(1,3,3);
    [nTraces3, nTrans3,tdpData] = tdplot_rtdata_v3('noPEP120.qub.dwt','noPEP120.traces',2,nTrans1);
    tplot( tdpData, tdpOptions{:} );
    xlim([-0.15 1.0]); xlabel('Initial FRET');
    ylim([-0.15 1.0]); ylabel('Final FRET');
    title('NoPEP Events');
    str1='T =  ';
    str2=num2str(nTrans3);
    str3='T/N =  ';
    str4=num2str(nTrans3/nTraces3);
    strA=strcat(str1,str2);
    strB=strcat(str3,str4);
    annotation('textbox',[0.825 0.730 0.074 0.05],'BackgroundColor',[1,1,1],...
    'string',strA,'FontSize',8);
    annotation('textbox',[0.825 0.675 0.074 0.05],'BackgroundColor',[1,1,1],...
    'string',strB,'FontSize',8);
end

%% 11. Make dwell time histogram.

disp('11. Making Dwell Time Histograms...');
figure;
set(gcf,'Position',[146 371 1037 404]);

if file_allMol_ac120>0
    subplot(1,3,1);
    dwtHst=popstate_itr('allMol_ac120.qub.dwt');
    plot(dwtHst(:,1),dwtHst(:,2),'ko',dwtHst(:,1),dwtHst(:,3),'go',dwtHst(:,1),dwtHst(:,4),'bo',dwtHst(:,1),dwtHst(:,5),'ro');
    xlim([-0.15 5.0]); xlabel('Time [s]');
    ylim([-0.05 1.05]); ylabel('Fraction');
    title('All Molecules');
end

if file_PEP120>0
    subplot(1,3,2);
    dwtHst=popstate_itr('PEP120.qub.dwt');
    plot(dwtHst(:,1),dwtHst(:,2),'ko',dwtHst(:,1),dwtHst(:,3),'go',dwtHst(:,1),dwtHst(:,4),'bo',dwtHst(:,1),dwtHst(:,5),'ro');
    xlim([-0.15 5.0]); xlabel('Time [s]');
    ylim([-0.05 1.05]); ylabel('Fraction');
    title('PEP Events');
end

if file_noPEP120>0
    subplot(1,3,3);
    dwtHst=popstate_itr('noPEP120.qub.dwt');
    plot(dwtHst(:,1),dwtHst(:,2),'ko',dwtHst(:,1),dwtHst(:,3),'go',dwtHst(:,1),dwtHst(:,4),'bo',dwtHst(:,1),dwtHst(:,5),'ro');
    xlim([-0.15 5.0]); xlabel('Time [s]');
    ylim([-0.05 1.05]); ylabel('Fraction');
    title('NoPep Events');
end
end



