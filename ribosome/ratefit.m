function results = ratefit(seg_fname, dwt_fname, bMakeGUI)
%RATEFIT  Derives average rates from filtered QuB MIL results
% 
%   RATES = RATEFIT(DIR, S, MODEL, SAMPLE_NAMES, STATE_NAMES)
%   
%   Returns a structure containing RATES from the QuB MIL (seperately)
%   results matrix 

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% TODO: split exponential fitting into new function: exp_ratefit
% Use input parser


% 4-state model (tRNA-tRNA FRET)
model = [ 4 3 ; 3 4 ; 3 2 ; 2 3 ; 4 2 ; 2 4 ];
state_names = {'PB','H2','H1','C'};

% 3-state model (tRNA-tRNA FRET)
% model = [ 2 3 ; 3 2 ];
% state_names = {'PB','H','C'};



%---- USER TUNABLE PARAMETERS ----

% Rate histogram binning parameters
BINSIZE=0.3;
logratesaxis=-3:BINSIZE:5;

% Histogram y-axis maximum
maxhist = 0.35;

% FILTERING PREFERENCES

% Minimum number of transitions to use rate in analysis
minTrans = 4;

% Minimum number of frames (trace length)
minLT = 120;

bMakeGUI = false;


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



% Prompt user for locations of files
if nargin<2,
    % rates file
    [datafile,datapath] = uigetfile( ...
                    '*.txt', 'Choose a QuB MIL-seperately results file:' );
    if datafile==0, return; end  %user hit 'cancel'
                
    seg_fname = [datapath filesep datafile];
    
    % idealization file
    [datafile,datapath] = uigetfile( ...
                    {'*.dwt'}, 'Choose a QuB Idealization file:' );
    disp(datafile);

    dwt_fname = [datapath filesep datafile];
    if datafile==0, return; end  %user hit 'cancel'
end

nSamples = 1;



%% Load segment files

% Load rates from QuB, filtering results
rates      = cell(nSamples,nRates);   % lists of MIL rates from QuB
rate_hist  = cell(nSamples,nRates);   % histogram of rates
rate_means = zeros(nSamples,nRates);  % means of guassian fit 
rate_stds  = zeros(nSamples,nRates);  % standard deviations
% special_rates = zeros(nSamples,4);    % locking,unlocking,kAA,kAP
fits = cell(nSamples,nRates);

for i=1:nSamples,

    assert( exist(seg_fname,'file')==2, ...
        sprintf('ERROR: Segments file %d does not exist',i) );
    
    % Read segments file, skipping first 2 columns and header row
    segments = dlmread(seg_fname,'\t',1,2);
    nTraces = size(segments,1);
    
    % read column header names
    fid  = fopen(seg_fname);
    cols = strread( fgetl(fid), '%s','delimiter','\t' );
    fclose(fid);
    
    % Get indeces of requested transitions into MIL results matrix
    rate_cols = zeros(1,size(model,1));
    
    for j=1:length(rate_cols)
        name = sprintf( 'k0 %d->%d', model(j,:) );
        rate_cols(j) = strmatch( name, cols );
    end
    
    r = log( segments( :, rate_cols ) );
    rates(i,:) = num2cell(r,1);
    
    
    % Count the number of transitions made (nSegments x nRates)
    assert( exist(dwt_fname,'file')==2, ...
        sprintf('ERROR: Idealization file %d does not exist',i) );
    
    ntrans = CountEvents( dwt_fname, model, size(segments,1) );
    [dwells,sampling] = loadDWT( dwt_fname );
    
    
    % Generate a hash table of the results matrix;
    % used for gathering filtering criteria (spc, ErrorCode)
    %h = cell2struct( num2cell(data), strrep(cols,' ','_'), 2);
    
    % Calculate length of each trace
    lt = zeros( nTraces,1 );
    
    for n=1:nTraces,
        lt(n) = sum( dwells{n}(:,2) );  %in frames
    end
    
    % CRITERIA FOR REMOVING MOLECULES
    % Remove rates for segments that did not make that transition.
    % Also remove segments that didn't converge.
    ncolErrorCode = strmatch( 'ErrorCode', cols );
    converged = segments(:,ncolErrorCode)==0;
    
    for j=1:nRates,
        data = rates{i,j};
    
        bad = ntrans(:,j)<minTrans | ~converged | isnan(data);
        %also use h.LL_per_event < -2
        
        rates{i,j} = data(~bad);
        disp( sum(~bad) );
    end

    disp( sprintf('N=%d', sum(converged)) );
    
    
    % Fit rate distributions to guassians
    for j = 1:nRates,
        
        data = rates{i,j};
        rate_hist{i,j} = hist(data,logratesaxis) / numel(data);
        fits{i,j} = fit( logratesaxis',rate_hist{i,j}', ...
                         'gauss1', 'Startpoint',[0.15,2,1.5] );

        %NOTE: this results in the MEDIAN of the distribution,
        %NOT THE MEAN! This will result in systematic overestimation
        %of kinetic parameters unless corrected...
        rate_means(i,j) = exp(fits{i,j}.b1);
        rate_stds(i,j)  = exp(fits{i,j}.c1);
    end
    
    % Find values for locking and unlocking rates
%     k_unlock = rate_means(i,1)+rate_means(i,5);
%     k_lock = klocking( rate_means(i,:) );
%     k_AA = rate_means(i,3) + rate_means(i,2);  %H1->H2 + H1->C
%     k_AP = kAP( rate_means(i,:) );
    
%     disp( sprintf( '%6s: %.2f (lock) %.2f (unlock) %.2f (kA/A) %.2f (kA/P) /s', ...
%           samples{i}, k_lock, k_unlock, k_AA, k_AP ) ...
%     ); 
%     special_rates(i,:) = [k_unlock k_lock k_AA k_AP];
    
end  % for each sample

disp(' ');




% RETURN VALUES FOR FUNCTION
results.rates = rates;
results.rate_hist = rate_hist;
results.rate_means = rate_means;
results.rate_stds  = rate_stds;
% results.special_rates = special_rates;

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
    elseif chosen_rate<nRates
        set(next, 'Enable', 'on');
    end

    ShowRates();
end  %function ShowPrev

function ShowNext(h, eventdata)
    chosen_rate = chosen_rate+1;
    if chosen_rate==nRates,
        set(next, 'Enable', 'off');
    elseif chosen_rate>1
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

    for i=1:nSamples,

        % Plot histogram and gaussian fit
        %subplot(s1,s2,i);  % top-down
        bar( logratesaxis,rate_hist{i,chosen_rate} );
        hold on, plot(fits{i,chosen_rate}), hold off;
        axis(histaxis);
        legend off;
        
        % Plot labels, etc
        xlabel(gca,  sprintf('log( k_{%s} )', t{chosen_rate} ), 'FontSize',14  );
        ylabel(gca, 'Percent in bin', 'FontSize',14 );
%         title( ['{' sample_names{i} '}'], 'FontSize',16, 'FontWeight','bold' );
       
        s = sprintf('%.2f (%.2f)', rate_means(i,chosen_rate), rate_stds(i,chosen_rate));
        
        text( histaxis(1)+0.5,0.9*histaxis(4), s, ...
              'FontSize',16,'HorizontalAlignment','left', ...
              'VerticalAlignment','middle');
        

    end %for each sample
    
end %function ShowRates


end %function ratefit (main)














