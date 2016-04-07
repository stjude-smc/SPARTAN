function varargout = occtime(varargin)
%OCCTIME  Average state occupancy over time
%
%   [OCC,TIME] = OCCTIME(FILES) is the occupancy of each FRET 
%   state over time in the given .dwt file, with states listed across columns  
%   and time across rows. The time axis is in seconds.
%   If files is a cell array, OCC is a cell array, one per file.
%
%   [...] = occupancyTimecourse() will prompt the user for files to load.
%
%   occupancyTimecourse(...) with no outputs displays the data in a new figure.
%
%   occupancyTimecourse(AX,...) plots in the the scalar axes AX.

%   Copyright 2016 Cornell University All Rights Reserved.



params.nFrames = 300;  %number of frames to show
% params.hideZeroState = false;



%% Process input arguments
narginchk(0,2);
nargoutchk(0,2);
[varargout{1:nargout}] = deal([]);

if nargin>0 && all(ishghandle(varargin{1},'axes'))
    ax = varargin{1};
    args = varargin(2:end);
else
    ax = [];
    args = varargin;
end

switch numel(args)
    case 0
        files = getFiles('*.dwt');
    case 1
        files = args{1};
    case 2
        [files,inputParams] = args{:};
        params = mergestruct(params, inputParams);
end
cellinput = iscell(files);
% files = findDwt(files);

nFiles = numel(files);
if nFiles==0,  return;  end

nFrames = params.nFrames;

if ~isempty(ax) && numel(ax)~=nFiles,
    error('Input axes must match number of files');
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
hasTarget = ~isempty(ax);
if ~hasTarget,
    hFig = figure;
end
titles = trimtitles(files);

% Format list of states (number, mean FRET value) for legend.
lgtxt = cell(nStates,1);
for s=1:nStates,
    lgtxt{s} = sprintf('State %d (%.2f)\t',s,fret(s));
end

for i=1:nFiles,
    % Clear axes for new plot.
    if ~hasTarget
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
end

linkaxes(ax,'xy');

end %function occupancyTimecourse


%% Save to file



% [dwt,sampling,~,fretModel] = loadDWT(dwtFile);
% nStates = size(fretModel,1);
% 
% occupancy = zeros(nStates,traceLength);
% time = 0:sampling:(traceLength-1)*sampling;
% 
% for i=1:length(dwt)
%     currentTrace = dwt{i};
%     currentPosition = 1;
%     %allStates = currentTrace(:,1);
%     for j=1:size(currentTrace,1)
%         currentState = currentTrace(j,1);
%         currentLength = currentTrace(j,2);
%         endPosition = currentPosition + currentLength - 1;
%         occupancy(currentState,currentPosition:endPosition) = occupancy(currentState,currentPosition:endPosition) + 1;
%         currentPosition = endPosition+1;
%     end
% end
% 
% for i=1:size(occupancy,2)
%     occupancy(:,i) = occupancy(:,i)./sum(occupancy(:,i));
% end
% 
% [path,name,~] = fileparts(dwtFile);
% outFile = fullfile(path, [name '_stateOcc.txt']);
% dlmwrite(outFile,vertcat(time,occupancy));
% 
% end
