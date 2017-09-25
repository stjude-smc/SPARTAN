function dwellplots(varargin)
%dwellplots  Display dwell-time histograms
% 
%   dwellplots(FILES,PARAMS) displays log-scale dwell-time histograms for each
%   .dwt file in the cell array FILES (must all have same time resolution).
%   PARAMS is a struct array with optional arguments.
%
%   dwellplots(FIG,...) specifies a target figure to draw the plot in.
%   
%   For internal use only. Use dwellhist() instead.
%
%   See also: dwellhist, lifetime_exp, loadDwelltimes, removeBlinks.

%   Copyright 2007-2017 Cornell University All Rights Reserved.

% FIXME: if targetting other axes, this will alter the menus...

narginchk(1,3)
nargoutchk(0,0);


%% Process input arguments calculate histograms
if ishandle(varargin{1}),
    hFig = varargin{1};
    varargin = varargin(2:end);
else
    hFig = figure;
end
[dwtfilename,params] = varargin{:};

if numel(dwtfilename)==0, return; end
names = trimtitles(dwtfilename);


% Calculate histograms
set(hFig,'pointer','watch'); drawnow;
[dwellaxis,histograms,fits] = dwellhist(dwtfilename,params);



%% Display survival plots, one state per panel.

% Find a good zoom axis range for viewing all of the histograms.
h = [histograms{:}];
ymax = 1.1*max(h(:));

if params.logX,
    xmax = dwellaxis(end);
else
    xmax = dwellaxis(  find( sum(h>0.01,2), 1, 'last' )  );
end

% If expected mean dwell times provided, show them as fit lines.
% Here, calculate normalization constants and change histogram line style.
if isfield(params,'model')
    lineStyle = 'b.';
else
    lineStyle = 'b-';
end

% Choose ordinate label based on normalization
if params.logX
    switch params.normalize
        case {'none','off'}
            ordinate = 'Counts';
        case 'state'
            ordinate = 'Counts (%)';
        case 'file'
            ordinate = 'Counts (% of file)';
        case 'time'
            ordinate = 'Counts s^{-1}';
        otherwise
            error('Invalid normalization setting');
    end
else
    ordinate = 'Dwell Survival (%)';
end


% Draw survival plots and fit lines, if model parameters given.
[~,nStates] = size(histograms);
ax = zeros(nStates,1);

for state=1:nStates,
    ax(state) = subplot( nStates, 1, state, 'Parent',hFig );

    if params.logX,
        semilogx( ax(state), dwellaxis, [histograms{:,state}], lineStyle, 'LineWidth',2 );
    else
        plot( ax(state), dwellaxis, [histograms{:,state}], lineStyle, 'LineWidth',2 );
    end
    
    ylabel(ax(state), ordinate);
    xlim( ax(state), [dwellaxis(1) xmax] );
    ylim( ax(state), [0 ymax] );
    title(ax(state), sprintf('State %d',state) );
end

% Draw fit lines
if isfield(params,'model')
    for state=1:nStates,
        hold( ax(state), 'on' );
        plot( ax(state), dwellaxis, fits(:,state), 'r-' );
        hold( ax(state), 'off' );
    end
end

xlabel(ax(end), 'Time (s)');
legend(ax(end), names);
linkaxes(ax,'xy');

set(hFig,'pointer','arrow'); drawnow;


%% Add menu items for adjusting settings and saving output to file
prompt = {'Remove blinks:', 'Log scale:', 'Log bin size:', 'Normalization:'};
fields = {'removeBlinks', 'logX', 'dx', 'normalize'};
types{4} = {'none','state','file','time'};
cb = @(~,~)settingdlg(params,fields,prompt,types,@dwellplots,{hFig,dwtfilename});
output = [to_col(dwellaxis) horzcat(histograms{:})];

defaultFigLayout( hFig, @(~,~)dwellplots(getFiles('*.dwt'),params), ...
                      @(~,~)dwellplots(hFig,getFiles('*.dwt'),params), ...
                      {@exportTxt,dwtfilename,output}, ...
       {'Change settings...',cb; ...
        %'Reset settings',@(~,~)dwellplots(hFig,dwtfilename) ...  %FIXME!
        'Copy output',{@clipboardmat,output}}  );


end %function dwellplots





%% ------ Save results to file for plotting in Origin
function exportTxt(~,~,files,output)
% Callback function for the "save histograms" button in the histogram figure.

names = trimtitles(files);
nFiles = numel(names);
nStates = (size(output,2)-1)/nFiles;

% Ask the user for an output filename.
[f,p] = uiputfile('*.txt','Select output filename',[mfilename '.txt']);
if f==0, return; end  %user hit cancel.
outFilename = fullfile(p,f);


% Output header lines
fid = fopen(outFilename,'w');
fprintf(fid,'Time (s)');

for state=1:nStates,
    for i=1:nFiles
        fprintf(fid,'\tState%d %s',state,names{i});
    end
end
fprintf(fid,'\n');
fclose(fid);


% Output histogram data
dlmwrite(outFilename, output, 'delimiter','\t', '-append');


end


