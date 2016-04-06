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

%   Copyright 2007-2016 Cornell University All Rights Reserved.

narginchk(1,3)
nargoutchk(0,0);


%% Process input arguments and create display window
if ishandle(varargin{1}),
    hFig = varargin{1};
    varargin = varargin(2:end);
else
    hFig = figure;
end

switch numel(varargin)
    case 1
        dwtfilename = varargin{1};
    case 2
        [dwtfilename,params] = varargin{:};
end

if numel(dwtfilename)==0, return; end
names = trimtitles(dwtfilename);


% Add menu items for adjusting settings and saving output to file.
hTxtMenu = findall(hFig, 'tag', 'figMenuGenerateCode');
set(hTxtMenu, 'Label','Export as .txt', 'Callback',@dwellplots_save);

prompt = {'Remove blinks:', 'Log scale:', 'Log bin size:', 'Normalization:'};
fields = {'removeBlinks', 'logX', 'dx', 'normalize'};
% types  = {'logical','logical','double',{'none','state','file','time'}};
cb = @(~,~,~)settingsDialog(params,fields,prompt,@dwellplots,hFig,dwtfilename);

hEditMenu = findall(hFig, 'tag', 'figMenuEdit');
delete(allchild(hEditMenu));
uimenu('Label','Display settings...', 'Parent',hEditMenu, 'Callback',cb);



%% Create histograms
set(hFig,'pointer','watch'); drawnow;
[dwellaxis,histograms] = dwellhist(dwtfilename,params);

handles.names = names;
handles.dwellhist = [to_col(dwellaxis) horzcat(histograms{:})];
guidata(hFig,handles);  %save for later calls to saveDwelltimes()



%% Display survival plots, one state per panel.

% Find a good zoom axis range for viewing all of the histograms.
% if params.logX,
    xmax = dwellaxis(end);
% else
%     xmax = 4* max(meanTime(:));
% end
h = [histograms{:}];
ymax = max(h(:));

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
end

% Draw survival plots
[~,nStates] = size(histograms);
ax = zeros(nStates,1);

for state=1:nStates,
    ax(state) = subplot( nStates, 1, state, 'Parent',hFig );

    if params.logX,
        semilogx( ax(state), dwellaxis, [histograms{:,state}], '-', 'LineWidth',2 );
        ylabel(ax(state), ordinate);
    else
        plot( ax(state), dwellaxis, [histograms{:,state}], '-', 'LineWidth',2 );
        ylabel(ax(state), 'Dwell Survival (%)');
    end
    if state==nStates,
        xlabel(ax(state), 'Time (s)');
    end
    
    xlim( ax(state), [dwellaxis(1) xmax] );
    ylim( ax(state), [0 ymax] );
    title(ax(state), sprintf('State %d',state) );
end

legend(ax(end), names);
linkaxes(ax,'xy');

set(hFig,'pointer','arrow'); drawnow;


end %function dwellplots





%% ------ Save results to file for plotting in Origin
function dwellplots_save(hObject,~,~)
% Callback function for the "save histograms" button in the histogram figure.

% Get histogram data saved in the figure object.
handles = guidata(hObject);
dwellhist = handles.dwellhist;
names = handles.names;

nFiles = numel(names);
nStates = (size(dwellhist,2)-1)/nFiles;

% Ask the user for an output filename.
[f,p] = uiputfile('*.txt','Select output filename','dwellhist.txt');
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
dlmwrite(outFilename, dwellhist, 'delimiter','\t', '-append');


end


