 function results = ratefit_exp(samples, model, state_names, bMakeGUI)
%RATEFIT_EXP  Estimates
% 
%   R = RATEFIT_EXP(D, S, b) ...
%   Calculates 

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% TODO: split exponential fitting into new function: exp_ratefit
% Use input parser


%---- USER TUNABLE PARAMETERS ----

% Rate histogram binning parameters
BINSIZE=0.4;
logratesaxis=-3:BINSIZE:5;

% Histogram y-axis maximum
maxhist = 0.25;

% Size of the histogram grid
s1 = 2;
s2 = 3;

% FILTERING PREFERENCES

% Minimum number of transitions to use rate in analysis
minTrans = 2;


%---------------------------------


if ~exist('bMakeGUI','var'), bMakeGUI=1; end  


% Generate transition rate names
nRates = size(model,1);

t = cell(nRates,1);

for i=1:nRates
    if ~exist('state_names','var')
        t{i} = sprintf( '{%d \\rightarrow %d}', model(i,:) );
    else
        t{i} = sprintf( '{%s \\rightarrow %s}', state_names{model(i,:)} );
    end
end %for each rate



% BASE FILENAMES of the samples, relative to the base path (d).
% Does not include extension!
nSamples = numel(samples);

% Display names for each sample (optional)
if ~exist('sample_names','var'),
    sample_names = strrep(samples,'_','\_');
end


% EXTENSIONS for files used:
dwt_ext  = '.qub.dwt';  % QuB Idealization file
seg_ext  = '.qub_seg.txt';  % QuB Segments output (whole matrix)
% sel_ext  = '.qub_sel.txt';  % QuB Selection list (not used)
% data_ext = '.qub.txt';  % forQuB raw data (not used)



%% Load segment files

% Load rates from QuB, filtering results
rates      = cell(nSamples,nRates);   % lists of MIL rates from QuB
rate_hist  = cell(nSamples,nRates);   % histogram of rates
rate_means = zeros(nSamples,nRates);  % means of guassian fit 
rate_stds  = zeros(nSamples,nRates);  % standard deviations
special_rates = zeros(nSamples,4);    % locking,unlocking,kAA,kAP
fits = cell(nSamples,nRates);

for i=1:nSamples,

    % Read segments file, skipping first 2 columns and header row
    segments = [basedir samples{i} seg_ext];
    assert( exist(segments,'file')==2, ...
        sprintf('ERROR: Segments file %d does not exist',i) );
    
    data = dlmread(segments,'\t',1,2);
    
    % read column header names
    fid = fopen(segments);
    cols = strsplit( fgetl(fid), '\t' );
    fclose(fid);
    
    assert( length(cols) == size(data,2) );
    
    % Generate a hash table of the data
    h = hashtable;
    
    for j=1:length(cols)
       h.( cols{j} ) = data(:,j);
    end
    
    
    % Gather rates together into a matrix
    rate_cols = zeros(1,size(model,1));
    
    for j=1:length(rate_cols)
       
        name = sprintf( 'k0 %d->%d', model(j,:) );
        rate_cols(j) = strmatch( name, cols );
        
    end
    
    r = log( data( :, rate_cols ) );
    rates(i,:) = num2cell(r,1);
    
    
    % Count the number of transitions made (nSegments x nRates)
    dwtfilename = [basedir samples{i} dwt_ext];
    assert( exist(dwtfilename,'file')==2, ...
        sprintf('ERROR: Idealization file %d does not exist',i) );
    
    ntrans = CountEvents( dwtfilename, model, size(data,1) );
    
    % CRITERIA FOR REMOVING MOLECULES
    % Remove rates for segments that did not make that transition.
    % Also remove segments that didn't converge.
    for r=1:nRates,
        data = rates{i,r};
    
        bad = ntrans(:,r)<minTrans | h.ErrorCode~=0 | isnan(data);
        
        rates{i,r} = data(~bad);
    end

    disp( sprintf('N=%d for %s', sum(h.ErrorCode==0), sample_names{i}) );
    
    
    % Fit rate distributions to guassians
    for r = 1:nRates,
        
        data = rates{i,r};
        rate_hist{i,r} = hist(data,logratesaxis) / numel(data);
        fits{i,r} = fit( logratesaxis',rate_hist{i,r}', ...
                         'gauss1', 'Startpoint',[0.15,2,1.5] );

        rate_means(i,r) = exp(fits{i,r}.b1);
        rate_stds(i,r)  = exp(fits{i,r}.c1);
    end
    
    % Find values for locking and unlocking rates
    k_unlock = rate_means(i,1)+rate_means(i,5);
    k_lock = klocking( rate_means(i,:) );
    k_AA = rate_means(i,3) + rate_means(i,2);  %H1->H2 + H1->C
    k_AP = kAP( rate_means(i,:) );
    
%     disp( sprintf( '%6s: %.2f (lock) %.2f (unlock) %.2f (kA/A) %.2f (kA/P) /s', ...
%           samples{i}, k_lock, k_unlock, k_AA, k_AP ) ...
%     ); 
    special_rates(i,:) = [k_unlock k_lock k_AA k_AP];
    
end  % for each sample

disp(' ');




% RETURN VALUES FOR FUNCTION
results.rates = rates;
results.rate_hist = rate_hist;
results.rate_means = rate_means;
results.rate_stds  = rate_stds;
results.special_rates = special_rates;

results.N = size(rates);
for i=1:nSamples,
    for j=1:nRates,
        results.N(i,j) = length(rates{i,j});
    end
end




%% Create GUI environment for cycling through rates

function ShowPrev(h, eventdata)
    chosen_rate = chosen_rate-1;
    if chosen_rate==1,
        set(prev, 'Enable', 'off');
    elseif chosen_rate==nRates-1
        set(next, 'Enable', 'on');
    end

    ShowRates();
end  %function ShowPrev

function ShowNext(h, eventdata)
    chosen_rate = chosen_rate+1;
    if chosen_rate==nRates,
        set(next, 'Enable', 'off');
    elseif chosen_rate==2
        set(prev, 'Enable', 'on');
    end

    ShowRates();
end  %function ShowNext


if bMakeGUI,
    f1 = figure();
    
    % These commands break matlab
    % set(f1,'DefaultAxesFontSize',12);
    % set(f1,'DefaulTtextFontSize',14);

    prev = uicontrol(f1, 'String', 'Prev', 'Callback', @ShowPrev, ...
        'Units', 'pixels', 'Position', [20 10 80 30]);
    next = uicontrol(f1, 'String', 'Next', 'Callback', @ShowNext, ...
        'Units', 'pixels', 'Position', [100 10 80 30]);

    chosen_rate = 1;
    set(prev, 'Enable', 'off');
    ShowRates  %display initial view
end %if bMakeGUI


%% 
function ShowRates()
    
    assert(chosen_rate>=1 & chosen_rate<=nRates, 'Invalid rate number');

    % Histogram bin parameters
    histaxis = [logratesaxis([1,end]) 0 maxhist]; %for plotting

    for i=1:numel(samples),

        % Plot histogram and gaussian fit
        subplot(s1,s2,i);  % top-down
        bar( logratesaxis,rate_hist{i,chosen_rate} );
        hold on, plot(fits{i,chosen_rate}), hold off;
        axis(histaxis);
        legend off;
        
        % Plot labels, etc
        xlabel(gca,  sprintf('log( k_{%s} )', t{chosen_rate} ), 'FontSize',14  );
        ylabel(gca, 'Percent in bin', 'FontSize',14 );
        title( ['{' sample_names{i} '}'], 'FontSize',16, 'FontWeight','bold' );
       
        s = sprintf('%.2f (%.2f)', rate_means(i,chosen_rate), rate_stds(i,chosen_rate));
        
        text( histaxis(1)+0.5,0.9*histaxis(4), s, ...
              'FontSize',16,'HorizontalAlignment','left', ...
              'VerticalAlignment','middle');
        

    end %for each sample
    
end %function ShowRates




%% Plot dwell time (state occupancy) histogram
function PlotDwellTime()
    % pop = zeros(numel(samples),4);  %4-state model, incl blinking

    dwellaxis = (1:1:100)-1;
    maxdwell = 1500;
    % histaxis = [1 35 0 0.3];

    % t = ['k12', 'k21', 'k23', 'k32', 'k13', 'k31'];
    before = [1 2 2 3 1 3];
    after  = [2 1 3 2 3 1];

    % figure();
    picked_sample = 3;
    i = picked_sample;
    disp( samples{i} );


    dwellhist = zeros( nRates, numel(dwellaxis) );

    % Load idealization file from QuB
    dwtfilename = [basedir samples{i} dwt_ext];
    dwells = LoadDWT( dwtfilename );

    for j=1:nRates, % for each transition in model

        trans_time = [];

        % Create a histogram of dwell times before each transition
        for seg=1:size(dwells,1),  % for each segment in dwelltime file
            [states,times] = dwells{seg,:};

            trans_time = [trans_time ; times( ...
                     states(1:end-1,1)==before(j) & ...
                     states(2:end,1)==after(j) ...
                     )];
        end

        data = histc(trans_time, 0:maxdwell);
        survival = sum(data) - cumsum(data); %num still in state after x frames
        dwellhist(j,:) = survival(dwellaxis+1);

        % Plot the exponential
        subplot(3,2,j);
        histplot = plot( dwellaxis, dwellhist(j,:), 'k.' );
%         set(gca,'XScale', 'log')
        xlim( [0 100] ); 
        title( t(j,:) );
    %     xlabel('Dwell time (frames)');
    %     ylabel('Percent in bin');
        hold on;

        % Fit the dwell times to an exponential (only amp varies)
        x = dwellaxis;
        y = dwellhist(j,:);

        [result1,goodness1] = fit( x', y', 'exp1' );
        exp1 = plot(result1, 'r');

        [result2,goodness2] = fit( x', y', 'exp2' );
        exp2 = plot(result2, 'b');

        % Also plot MIL rates
        MIL_y = y(1)* exp( -rate_means(i,j)*x*0.025 );  %last factor: ms to frames
        mil = plot( x,MIL_y, 'g' );

        legend( [exp1,exp2,mil], 'Single', 'Double', 'MIL' );
        hold off;

        % Display all rates and errors
        disp( sprintf('%s   single=%5.2f  MIL=%5.2f   dbl=%5.2f %5.2f (%4.2f %4.2f)', ...
              t(j,:), -result1.b/0.025, rate_means(i,j), ...
              -result2.b/0.025, -result2.d/0.025, ...
              result2.a, result2.c ) );

    end %for each rate

    disp(' ');
    
end %function PlotDwellTime




end %function ratefit (main)














