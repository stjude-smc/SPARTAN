function varargout = occtime(varargin)
%OCCTIME  Average state occupancy over time
%
%   [OCC,TIME] = OCCTIME(FILES) is the occupancy of each FRET 
%   state over time in the given .dwt file, with states listed across columns  
%   and time across rows. The time axis is in seconds.
%   If files is a cell array, OCC is a cell array, one per file.
%
%   [...] = occtime() will prompt the user for files to load.
%
%   occtime(...) with no outputs displays the data in a new figure.
%
%   occtime(AX,...) plots in the the axes AX (one per file)

%   Copyright 2016 Cornell University All Rights Reserved.


persistent params;
if isempty(params),
    params.nFrames = 300;  %number of frames to show
    params.hideZeroState = true;
end


%% Process input arguments
narginchk(0,3);
nargoutchk(0,2);
[varargout{1:nargout}] = deal([]);
[ax,hFig] = deal([]);

if nargin>0,
    if all(ishghandle(varargin{1},'axes'))
        ax = varargin{1};
        hFig = get(ax(1),'Parent');
        varargin = varargin(2:end);
    elseif all(ishghandle(varargin{1},'figure'))
        hFig = varargin{1};
        varargin = varargin(2:end);
    end
end

switch numel(varargin)
    case 0
        files = getFiles('*.dwt');
    case 1
        files = varargin{1};
    case 2
        [files,inputParams] = varargin{:};
        params = mergestruct(params, inputParams);
end
cellinput = iscell(files);
if ~iscell(files), files={files}; end
files = findDwt(files,'raiseError');

nFiles = numel(files);
if nFiles==0,  return;  end

nFrames = params.nFrames;

if ~isempty(ax) && numel(ax)~=nFiles,
    error('Input axes must match number of files');
end

% Create figure, 
if isempty(hFig) && nargout==0
    hFig = figure;
end
if ~isempty(hFig),
    set(hFig,'pointer','watch'); drawnow;
end



%% Calculate occupancy
occupancy = cell(nFiles,1);

for f=1:nFiles,
    % Load idealization
    [dwt,sampling,~,fretModel] = loadDWT(files{f});
    idl = dwtToIdl(dwt);
    idl = idl(:,1:nFrames);

    if f==1,
        time = (0:(nFrames-1))' *sampling/1000; %seconds
        fret = fretModel(:,1);
        nStates = numel(fret);
    else
        %FIXME check all files have the same sampling and nStates.
    end
    
    % Calculate occupancy, summing all rows to count traces at each timepoint.
    occ = zeros(nFrames,nStates);
    for s=1:nStates,
        occ(:,s) = sum(idl==s)';
    end
    
    if params.hideZeroState
        occ = occ(:,2:end);
    end
    
    occupancy{f} = 100*bsxfun(@rdivide,occ,sum(occ,2));
end

% Return in the same format is input (matrix in, matrix out).
if nargout>0
    if ~cellinput,
        output = {occupancy{1},time};
    else
        output = {occupancy,time};
    end
    [varargout{1:nargout}] = output{1:nargout};
    if isempty(ax), return; end
end



%% Display plots with occupancy
titles = trimtitles(files);

h = [occupancy{:}];
ymax = 1.1*max(h(:));

% Format list of states (number, mean FRET value) for legend.
lgtxt = cell(nStates,1);
for s=1:nStates,
    lgtxt{s} = sprintf('%d (%.2f)\t',s,fret(s));
end
if params.hideZeroState
    lgtxt = lgtxt(2:end);
end

doCreatePlot = isempty(ax);

for i=1:nFiles,
    % Clear axes for new plot.
    if doCreatePlot
        ax(i) = subplot(1,nFiles,i, 'Parent',hFig);
    else
        newplot(ax(i));
    end
    
    % Plot occupancy
    plot(ax(i), time, occupancy{i});
    
    if i==1,
        xlabel(ax(i), 'Time (s)');
        ylabel(ax(i), 'Occupancy (%)');
        legend(ax(i), lgtxt);
    end

    title(ax(i), titles{i});
    xlim(ax(i), [time(1) time(end)]);
    ylim(ax(i), [0 ymax]);
end

linkaxes(ax,'xy');
set(hFig,'pointer','arrow'); drawnow;



%% Add menus to change settings, get data, open new plots, etc.
txtout = [time occupancy{:}];

defaultFigLayout( hFig, @(~,~)occtime(getFiles('*.dwt'),params), ...  %File->New callback
                   @(~,~)occtime(hFig,getFiles('*.dwt'),params), ...  %File->Open callaback
                   {@exportTxt,files,txtout}, ...                     %Export callaback
      { 'Change settings...',@(~,~)settingdlg(params,@occtime,{ax,files}); ...
        'Reset settings',    @(~,~)occtime(ax,files); ...
        'Copy values',      {@clipboardmat,txtout}  }  );


end %function occupancyTimecourse


%  =========================================================================  %
%% Save the output to file.
function exportTxt(~,~,files,output)

nFiles = numel(files);
nStates = (size(output,2)-1)/nFiles;

if numel(files)==1,
    [p,f] = fileparts(files{1});
    outputFilename = fullfile(p, [f '_' mfilename '.txt']);
else
    outputFilename = [mfilename '.txt'];
end

[f,p] = uiputfile('*.txt', [mfilename ': save output'], outputFilename);
if f==0, return; end  %user hit cancel.
outFilename = fullfile(p,f);


% Output header lines
fid = fopen(outFilename,'w');
fprintf(fid,'Time (s)');
names = trimtitles(files);

for i=1:nFiles
    for state=1:nStates,
        fprintf(fid,'\t%s (State %d)',names{i},state);
    end
end
fprintf(fid,'\n');
fclose(fid);


% Output histogram data
dlmwrite(outFilename, output, 'delimiter','\t', '-append');


end %function exportTxt


