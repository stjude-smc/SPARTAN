function [lifetimes,fits,totalTimes,dwellaxis,dwellhist] = lifetime_exp( dwtfilename, inputParams )
%LIFETIME_EXP  Estimates median lifetime in each state
% 
%   R = LIFETIME_EXP(FILE)
%   Returns the lifetimes of each state in ascending order using the data
%   in FILE (filename of a .dwt file from QuB).  If FILE is a cell array
%   with multiple entries, lifetimes for all files will be returned as a
%   matrix, with files in rows and states in columns.
%   


%---- USER TUNABLE PARAMETERS ----

params.bMakeGUI = 1;
params.useCorrectedDwelltimes = 1;  % merge blinks into previous dwell

% Option to remove dwells whose durations are unknown because they are
% cropped by the start of measurement, blinking, and photobleaching, resp.
params.dropFirstDwell = 0;
params.dropLastDwell  = 0;  %take care of in tIdealize now...
params.dropDarkDwells = 0;  %before and after dwell in dark state

% colors for statehist, in order
params.colors = [ 0 0 0   ; ... % black
               0        0.7500   0.7500 ; ... % cyan
               0.7500   0        0.7500 ; ... % purple
               1 0 0   ; ... % red
               0.7500   0.7500   0 ; ...    % yellow
               0 0.5 0 ]; % green

params.fitSingle = true;  %otherwise, double exponential fitting...

params.plotFits = 1;

% Merge options, giving the user's options precedence.
if nargin>1,
    params = catstruct( params, inputParams );
end

%---------------------------------

% if no files give, prompt the use
if ~exist('dwtfilename','var'),
    disp('Select DWT files, hit cancel when finished');
    dwtfilename = getFiles('*.dwt','Choose a DWT file:');
end

if size(dwtfilename,1) == 0, return; end

% if just one filename given, convert to cell array
if ~iscell( dwtfilename ),
    dwtfilename = {dwtfilename};
end

nFiles = numel(dwtfilename);



%% Load segment files

% Get number of states..
[dwells,sampling,offsets,model] = LoadDWT( dwtfilename{1} );
nStates = numel(model)/2;
clear dwells; clear offsets;

%
lifetimes = zeros(nFiles,nStates);   % average lifetimes (fit)
totalTimes= zeros(nFiles,nStates);
dwellhist = cell(nFiles,nStates);   % histogram of dwell times
fits      = cell(nFiles,nStates);    % fit structures for plotting

for i=1:nFiles,
    
    assert( exist(dwtfilename{i},'file')==2, ...
            sprintf('ERROR: No such file: %s',dwtfilename{i}) );
    
    % Load DWT file (states in columns)
    if params.useCorrectedDwelltimes
        [dwells,sampling,model] = correctedDwelltimes( dwtfilename{i} );
    else
        [dwells,sampling,model] = loadDwelltimes( dwtfilename{i} );
    end
    assert( numel(dwells) == nStates );
        
    % Create a survival plot for each state
    dwellaxis = (0:1:5000)*sampling /1000;
    %dwellhist = zeros( nStates,numel(dwellaxis) );
    
    for j=2:nStates,
        if isempty( dwells{j} ),
            continue;
        end
        
        % Create survival plot
        times = dwells{j};
        totalTimes(i,j) = sum(times);
        data = histc( times./1000, dwellaxis );
        
        survival = sum(data) - cumsum(data);
        survival = survival/survival(1);
        dwellhist{i,j} = survival;
        
        % Fit survival plot to a single exponential
        % (truncated 0-tail to avoid fitting errors)
        plen = min( numel(dwellaxis), max(times)/sampling );
        x = dwellaxis(1:plen);
        y = survival(1:plen);
        
        if params.fitSingle,
            if j>1,
                result1 = fit( x', y, 'exp1', 'StartPoint',[1 -1/mean(times)/1000] );
            else
                result1 = fit( x', y, 'exp1' );
            end
            
            coefs = coeffvalues(result1);
            if j>1 && coefs(1)<0.7,
               warning(['Not a great exponential fit: ' dwtfilename{i} ]);
            end
        else
            result1 = fit( x', y, 'exp2' );
        end
        
        % Find weighted average of time constants.
        % If a single exponential, just retrieves the time constant.
        fits{i,j} = result1;
        coefs = coeffvalues(result1);
        weights = coefs(1:2:end);
        weights = weights./sum(weights);

        % Record mean lifetime
        weightedAvg = mean( coefs(2:2:end).*weights );
        lifetimes(i,j) = -1/weightedAvg;
    end

end  % for each sample


if nFiles>size(params.colors,1),
    params.colors = zeros(nFiles,3);
end

% Make a seperate figure that combines the plots for each inspection.
if params.plotFits,
    h1 = figure();

    set(h1,'DefaultAxesColorOrder',params.colors);

    nrows = nStates-1;
    ncols = nFiles+1;
    
    for i=1:nFiles,

        % Plot distribution and fit of each state
        for j=2:nStates
            % Plot data and fit
            subplot( nrows,ncols, (ncols*(j-2))+i );
            cla;
            plot( dwellaxis, dwellhist{i,j}, 'k.' ); hold on;
            hp = plot( fits{i,j}, 'r-');
            set(hp,'LineWidth',2);
            xlim( [0 4*mean(lifetimes(:,2))] );
            ylim( [0 1] );
            legend off;

            if j==2,
                xlabel(gca, 'Time (sec)', 'FontSize',14  );
            end
            if i==1,
                ylabel(gca, 'Dwell Count (%)', 'FontSize',14 );
            end

            % Plot this also in the overlay plot at end
    %         subplot( nrows,ncols, ncols*(j-1) );
    %         
    %         plot( dwellaxis, dwellhist{i,j}, '.', 'MarkerEdgeColor',params.colors(i,:) );
    %         xlim( [0 1.5] );
    %         ylim( [0 1] );
    %         hold on;
    % 
    %         if j==nStates,
    %             xlabel(gca, 'Time (sec)', 'FontSize',14);
    %         end
        end

    end
end

h2 = figure;
set(h2,'defaultaxesfontsize',16);
set(h2,'defaulttextfontsize',18);


% Make figure that overlays decay curves for direct comparison
ax = [];

output = dwellaxis';

if ~isfield(params,'colors'),
    params.colors = colormap;
    nlevels = size(params.colors,1);
    params.colors = colors(1:round(nlevels/nFiles):end, :);
end

for i=1:nFiles,

    for j=2:nStates

        % Plot this also in the overlay plot at end
        ax(j-1) = subplot( nStates-1,1, j-1 );
        
        plot( dwellaxis, dwellhist{i,j}, '-', 'Color',params.colors(i,:), 'LineWidth',2 );
%         plot( dwellaxis, dwellhist{i,j}, 'r-', 'LineWidth',2 );
        xlim( [0 60] );
        ylim( [0 1] );
        set(gca,'ytick',[0:10:60]);
        hold on;
        
        title( ['State ' num2str(j) ' Lifetime'] );
        
%         if j==2,
            ylabel(gca, 'Dwell Count (%)' );
%         else
%             set(gca,'yticklabel',[]);
%         end
        if j==nStates,
            xlabel(gca, 'Time (sec)');
        end
    end
end

legend;
% linkaxes(ax,'y');


%------ Save results to file for plotting in Origin

% Output header lines
fid = fopen('lifetime_exp.txt','w');
fprintf(fid,'Time (sec)\t');

for j=2:nStates,
    for i=1:nFiles
        [p,name] = fileparts( dwtfilename{i} );
        fprintf(fid,'State%d %s\t',j,name);
        output = [output dwellhist{i,j}];
    end
end
fprintf(fid,'\n');

% Output data
for i=1:size(output,1),
    fprintf(fid,'%d\t',output(i,:));
    fprintf(fid,'\n');
end

fclose(fid);


%% Create GUI environment for cycling through rates

function ShowPrev(h, eventdata)
    chosenFile = chosenFile-1;
    if chosenFile==1,
        set(prev, 'Enable', 'off');
    elseif chosenFile==nFiles-1
        set(next, 'Enable', 'on');
    end

    ShowRates();
end  %function ShowPrev

function ShowNext(h, eventdata)
    chosenFile = chosenFile+1;
    if chosenFile==nFiles,
        set(next, 'Enable', 'off');
    elseif chosenFile==2
        set(prev, 'Enable', 'on');
    end

    ShowRates();
end  %function ShowNext


if false && bMakeGUI,
    f1 = figure();
    
    % These commands break matlab
    % set(f1,'DefaultAxesFontSize',12);
    % set(f1,'DefaulTtextFontSize',14);

    prev = uicontrol(f1, 'String', 'Prev', 'Callback', @ShowPrev, ...
        'Units', 'pixels', 'Position', [20 10 80 30]);
    next = uicontrol(f1, 'String', 'Next', 'Callback', @ShowNext, ...
        'Units', 'pixels', 'Position', [100 10 80 30]);

    chosenFile = 1;
    set(prev, 'Enable', 'off');
    ShowRates  %display initial view
end %if bMakeGUI


function ShowRates()
    
    assert(chosenFile>=1 & chosenFile<=nFiles, 'Invalid file number');

    % Histogram bin parameters
    %histaxis = [logratesaxis([1,end]) 0 maxhist]; %for plotting

    for j=2:nStates,
        
        if isempty( dwellhist{chosenFile,j} )
            continue;
        end
        
        % Plot data and fit
        subplot(1,nStates-1,j-1);
        cla;
        plot( dwellaxis, dwellhist{chosenFile,j}, 'k.' ); hold on;
        hp = plot( fits{chosenFile,j}, 'r');
        set(hp,'LineWidth',2);
        %axis( dwellaxis );
        xlim( [0 12] );
        legend off;
        
        % Plot labels, etc
        xlabel(gca, 'Time (ms)', 'FontSize',14  );
        ylabel(gca, 'Num. Still Alive', 'FontSize',14 );
       
%         text( histaxis(1)+0.5,0.9*histaxis(4), s, ...
%               'FontSize',16,'HorizontalAlignment','left', ...
%               'VerticalAlignment','middle');

    end %for each sample
    
    if exist('titles','var'),
        title( titles{chosenFile} );
    end
    
end %function ShowRates






end %function...




