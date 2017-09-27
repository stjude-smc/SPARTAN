function [lifetimes,rsquare,fits] = lifetime_exp( dwtfilename, inputParams )
% dwellhist  Survival plots of state dwell-times
% 
%   [LIFETIMES,RSQUARE,FITS] = lifetime_exp(FILES) creates survival plots of
%   state dwell-times for each state in each .dwt file in the cell array FILES.
%   These are then fit to exponential decay functions to estimate the decay
%   time constants LIFETIMES (files in rows, states in columns). The adjusted
%   R-squared goodness of fit scores RSQUARE and the fit objects FITS are also
%   returned.
%
%   [...] = lifetime_exp(...,PARAMS) specifies optional parameters in the 
%   struct PARAMS to control how the histograms are made (true/false values).
%   Most are used for dwellhist(), which makes the survival plots.
%
%      'removeBlinks': Ignore dwells in the dark state, assumed to be state 1.
%                      Dwells broken up by such blinks are merged, with the
%                      time during the blink added to surrounding dwells.
%
%      'fitSingle':    Fit a single (true) or double (false) exponential.
% 
%      'plotFits':     Display all of the individual fits.
%
%   See also: dwellhist, loadDwelltimes, removeBlinks.

%   Copyright 2007-2015 Cornell University All Rights Reserved.



%% ---- USER TUNABLE PARAMETERS ----
params.removeBlinks = true;  % merge blinks into previous dwell

params.logX = false;  %no Sine/Sigworth transform.
params.fitSingle = true;  %otherwise, double exponential fitting...

params.plotFits = true;

% Merge options, giving the user's options precedence.
if nargin>1,
    params = mergestruct( params, inputParams );
end



%% Prompt user for file names if not given.
if nargin<1,
    disp('Select DWT files, hit cancel when finished');
    dwtfilename = getFiles('*.dwt','Choose dwell-time files');
end

if ischar(dwtfilename),
    dwtfilename = {dwtfilename};
end

nFiles = numel(dwtfilename);
if nFiles==0, return; end
names = trimtitles(dwtfilename);



%% Load dwell-times and calculate survival plots
[dwellaxis,survival,meanTimes] = dwellhist(dwtfilename, params);
nStates = size(survival,2);



%% Fit survival plots to exponential functions
lifetimes = zeros(nFiles,nStates);   % average lifetimes (fit)
rsquare = zeros(nFiles,nStates);    % R-squared values for each fit.
fits = cell(nFiles,nStates);    % fit structures for plotting

for file=1:nFiles,
    for state=1:nStates,
        % Truncate the survival plot to de-emphasize the tail for fitting.
        y = survival{file,state};
        plen = find(y>0.01,1,'last');
        x = dwellaxis(1:plen)';
        y = y(1:plen);
        
        % Fit the survival plot to an exponetial
        if params.fitSingle,
            [fits{file,state},gof] = fit( x, y, 'exp1', ...
                                 'StartPoint', [1 -1/meanTimes(file,state)] );
        else
            [fits{file,state},gof] = fit( x, y, 'exp2' );
        end
        
        % Find weighted average of time constants.
        % If a single exponential, just retrieves the time constant.
        coefs = coeffvalues( fits{file,state} );
        weights = coefs(1:2:end);
        weights = weights./sum(weights);

        % Record mean lifetime
        weightedAvg = mean( coefs(2:2:end).*weights );
        lifetimes(file,state) = -1/weightedAvg;
        rsquare(file,state) = gof.adjrsquare;
    end

end  % for each sample

if ~params.plotFits, return; end



%% Display fits

% Save the histogram data in the figure and add a button so the data
% can be saved to file by the user.
hFig = figure;
output = [to_col(dwellaxis) horzcat(survival{:})];
setappdata(hFig,'survival',output);
setappdata(hFig,'names',names);

nrows = nStates;
ncols = nFiles+1;
xend = max( 3.5*lifetimes(:) );

ax = zeros(nrows,ncols);

for i=1:nFiles+1,    
    for j=1:nStates,
        
        % Plot survival plot fits.
        if i<=nFiles,
            %idxCol = nFiles*(j-1)+i;
            %y = survival(:,idxCol);
            y = survival{i,j};

            % Plot data and fit
            ax(i,j) = subplot( nrows,ncols, (ncols*(j-1))+i, 'Parent',hFig );
            plot( ax(i,j), dwellaxis, y, 'k.' ); hold on;
            hp = plot( fits{i,j}, 'r-');  %FIXME: ax(i,j)
            set(hp,'LineWidth',2);
            legend(ax(i,j), 'off');
            
            if j==1, %first row
                title(names{i});
            end
        
        % Overlay all survival plots for direct comparison
        else
            %idxCol = nFiles*(j-1)+(1:nFiles);
            %y = survival(:,idxCol);
            y = [survival{:,j}];

            % Plot data and fit
            ax(i,j) = subplot( nrows,ncols, (ncols*(j-1))+i, 'Parent',hFig );
            plot( ax(i,j), dwellaxis, y, 'LineWidth',2 );
            
            if j==1,
                legend(ax(i,j), names);
                title(ax(i,j), 'Overlay');
            end
        end
            
        % Set axis and plot titles, etc.
        xlim( ax(i,j), [0 xend] );
        ylim( ax(i,j), [0 1] );
        
        if j==nStates,
            xlabel(ax(i,j), 'Time (s)' );
        else
            xlabel(ax(i,j), '');
        end

        if i==1, %first column)
            ylabel( ax(i,j), sprintf('State %d',j) );
        else
            ylabel(ax(i,j), '');
        end
        
    end %for each state
end %for each file

linkaxes(ax,'xy');

% Add a control at the bottom of the GUI for saving the histograms to file.
uicontrol( 'Style','pushbutton', 'String','Save...', ...
           'Position',[15 15 75 30], 'Callback',@saveDwelltimes, ...
           'Parent',hFig );



end %function...



%% ------ Save results to file for plotting in Origin
function saveDwelltimes(hObject,~,~)
% Callback function for the "save histograms" button in the histogram figure.

% Get histogram data saved in the figure object.
hFig = get(hObject,'Parent');
dwellhist = getappdata(hFig,'survival');
names = getappdata(hFig,'names');

nFiles = numel(names);
nStates = (size(dwellhist,2)-1)/nFiles;


% Ask the user for an output filename.
[f,p] = uiputfile('*.txt','Select output filename',[mfilename '.txt']);
if f==0, return; end
outFilename = fullfile(p,f);


% Output header lines
fid = fopen(outFilename,'w');
fprintf(fid,'Time (s)\t');

for state=1:nStates,
    for i=1:nFiles
        fprintf(fid,'State%d %s\t',state,names{i});
    end
end
fprintf(fid,'\n');

% Output histogram data
for i=1:size(dwellhist,1),
    fprintf( fid, '%d\t', dwellhist(i,:) );
    fprintf( fid, '\n' );
end

fclose(fid);

end

