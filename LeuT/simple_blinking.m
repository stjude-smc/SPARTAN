%% Simulation of the effects of blinking on
%  LeuT kinetic parameter estimates
% 10/31/08
% Daniel Terry

traceID  = 1;   %keep track of the currently selected trace


%% --- 1. SIMULATE 40ms traces, long lived low FRET, short lived high

nTraces = 200;
framerate = 25;  %40ms
seed = 1272729;
traceLen = 1000;


%--- Structural kinetic model

% Average state lifetimes
tauO = 2;
tauC = 0.2;

% Emission probability parameters
%         open        closed
model = [ 0.52 0.07 ; 0.75 0.07 ];
mu = model(:,1);
sigma = model(:,2);
nStates = length(mu);

% Steady-state probabilities
pO = tauO/(tauC+tauO)
pC = tauC/(tauC+tauO)

% Rate matrix (simple 2-state model)
%     ->O   ->C
Qs = [ 0       1/tauO ;
      1/tauC  0        ];

p0s = [pO pC];



% --- Simulate the structural model to get raw data

[rawData,rawIdl] =  simulate( Qs, p0s, 0, model, nTraces, traceLen, ...
                                 framerate, seed );
% save( 'simple1.mat', 'rawData', 'rawIdl' );



%% --- Photophysical kinetic model

% Average state lifetimes
tauTOn =  30;   % average time before photobleaching
tauOn  =  1.25;   % average time before blink
tauOff =  0.25; % average time in the "blinked" state

% Emission probability parameters
stdDark = 0.05;
modelDark = [ 0.01 stdDark ; 0.01 stdDark ];

% Steady-state probabilities
pOn  = tauOn/(tauOn+tauOff);
pOff = tauOff/(tauOn+tauOff);

% Rate matrix (simple 2-state model)
%     ->off   ->on
Qp = [ 0         1/tauOff ;
       1/tauOn  0        ];
   
p0p = [0 1];
   
[blinkingData,blinkingIdl] =  simulate( Qp, p0p, 1/tauTOn, modelDark, ...
                                 nTraces, traceLen, framerate, (seed+7)/2 );

% There is a better way to do this
filteredData = rawData;
bi = blinkingIdl( :, 1:40:end );
filteredData(bi==1) = blinkingData(bi==1);

% saveTraces( 'noiseless_blinking2.qub.txt', 'qub', filteredData );


%% Plot traces for examination

i = 26;

% subplot( 2,1,1);
x  = (1:traceLen)*0.04;
xf = (1:size(rawIdl,2))*0.001;
% plot( xf, mu( rawIdl(i,:) ), 'r-', x, rawData(i,:), 'b-' );
plot( xf, mu( rawIdl(i,:) ), 'r-', x, filteredData(i,:), 'b-' );
xlim( [1 20] );
ylim( [-0.1 1] );

% subplot( 2,1,2);
% hist( filteredData(j,ll),0:0.02:0.8 )
% xlim( [0 0.8] );
% pbaspect( [1 1 1] );










%% INITIAL CONDITIONS FOR ALL PARAMETER ESTIMATION RUNS

nStates = 3;
p0_start = ones(1,nStates)/nStates;     % state initial probabilities

A_start       = ones(nStates)/nStates;       % transition probability matrix
x =2/framerate; y=1-((nStates-1)*x);              % initial rate=2/sec
A_start(:) = x;  A_start( find(eye(nStates)) ) = y;

mu_start = [0.01; mu]';
sigma_start = [0.05; sigma]';




% --- 2. PARAMETER ESTIMATES FOR Baum-Welch TOGETHER
% profile clear;

nRates   = nStates^2 - nStates;

% Get parameter estimates many times
% LL = cell(nRuns-1,1);
% A  = cell(nRuns-1,1);
% p0 = zeros(nRuns-1,nStates);
% ps = zeros(nRuns-1,nStates);

[LL,A,muE,sigmaE,p0,ps] = BWoptimize( ...
            filteredData, A_start, mu_start, sigma_start, ...
            p0_start,'FixMu', ones(1,nStates), 'FixSigma',ones(1,nStates) );


% 1->[2,3,4], 2->[1,3,4], 3->[1,2,4], 4->[1,2,3]
Qt = framerate.*A';
disp( [ Qt(2,1) Qt(1,2)+Qt(1,3) Qt(3,2) Qt(2,3) ] );
 
%%
Q_parts(i,:)  = Qt( [2 3 4   5 7 8   9 10 12   13 14 15] );
    

% Plot the observed distribution, noting actual value
% This will all be in a seperate function:
% [fits,fig] = plotSimResults( A,ps, A_true?, ps_true? )
% set(fig,'Name','Baum-Welch together');
% 
% subplot
% kCH1, kH1H2, kCH2
% kH1C, kH2H1, kH2C
% ps1,  ps2,   ps3

%% Fit each to get stderr, etc for comparison
% and mark the true mean
rates = trnaRates( Q_parts );
sprintf('%.2f  ', mean(rates) )
sprintf('%.2f  ', std(rates)  )

ps2 = ps(:,2:end);
ps2n = ps2./repmat( sum(ps2,2),1,3 );
sprintf('%.2f  ', mean(ps2n) )
sprintf('%.2f  ', std(ps2n)  )


% Save the results to file for records
save( 'bw_together2.mat', 'LL', 'A', 'p0', 'ps' );



%% Idealize using the BW estimations to understand why the
% estimates are so bad
profile clear;

% Use averages of the above parameter estimates for idealization
p0_fit = mean(p0);

A_fit  = zeros(nStates,nStates);
for i=1:length(A),
    A_fit = A_fit + A{i};
end
A_fit = A_fit./length(A);



% Find the viterbi sequence of each trace
[NT,len] = size(filteredData);

bwIdl = zeros( NT,len );

h = waitbar(0,'Idealizing traces...');
for i=1:NT
    obs = filteredData(i,:);
    Bx  = constructB(obs,mu,sigma);
    [LL, vPath, vLL] = forward_viterbi(p0_fit, A_fit, Bx');
    bwIdl(i,:) = vPath;
    
    waitbar(i/NT,h);
end
close(h);



%% Plot traces for examination

i = 94;

idxs = find(idx);
j = idxs(i);

ll  = 1:1000 /2; %show half the trace
ll2 = 1:25000;

subplot( 2,1,1);
stairs( ll2/25, mu( idl(j,ll2) ), 'k-' ); hold on;
stairs( ll, mu( bwIdl(i,ll) ), 'r-' );
plot( ll, rawdata(j,ll), 'b-' ); hold off;
xlim( [1 ll(end)] );
ylim( [-0.1 0.8] );

subplot( 2,1,2);
hist( filteredData(j,ll),0:0.02:0.8 )
xlim( [0 0.8] );
pbaspect( [1 1 1] );



%% --- 3. PARAMETER ESTIMATES FOR Baum-Welch 5 GROUPS



%% --- 4. PARAMETER ESTIMATES FOR Baum-Welch INVIVIDUALLY,FILTERED



%% --- 5. MIL TOGETHER

% Split data into 5 subsets and save them for QuB
for i=1:(nRuns-1)
    s = nPerRun*(i-1) +1;
    e = nPerRun*i;
    
    disp( [s e] );
    
    fname = sprintf('sim1_filtered_%d.qub.txt',i);
    
    saveTraces( fname, 'qub', filteredData(s:e,:) );
end


%% --- 6. MIL Seperately

dwt_fnames = { ...
    'sim1_filtered_1.qub.dwt', ...
    'sim1_filtered_2.qub.dwt', ...
    'sim1_filtered_3.qub.dwt', ...
    'sim1_filtered_4.qub.dwt', ...
    'sim1_filtered_5.qub.dwt'};

seg_fnames = { ...
    'sim1_filtered_1.qub_seg.txt', ...
    'sim1_filtered_2.qub_seg.txt', ...
    'sim1_filtered_3.qub_seg.txt', ...
    'sim1_filtered_4.qub_seg.txt', ...
    'sim1_filtered_5.qub_seg.txt'};

results = cell(0,1);

for i=1:length(dwt_fnames);
    results{i} = ratefit( seg_fnames{i}, dwt_fnames{i} );
end

%% Process the result

% Collect the results into a matrix
mils_mu    = zeros(length(results),6);
% mils_sigma = zeros(length(results),6);

for i=1:length(results)
    mils_mu(i,:)    = results{i}.rate_means;
%     mils_sigma(i,:) = results{i}.rate_stds;
end

% Find the mean of the results and stdev.
% rates = trnaRates( Q_parts );
sprintf('%.2f  ', mean(mils_mu) )
sprintf('%.2f  ', std(mils_mu)  )

% ps2 = ps(:,2:end);
% ps2n = ps2./repmat( sum(ps2,2),1,3 );
% sprintf('%.2f  ', mean(ps2n) )
% sprintf('%.2f  ', std(ps2n)  )






