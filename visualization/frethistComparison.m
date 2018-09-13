function varargout = frethistComparison(varargin)
%fresthistComparison  1D FRET histogram
%
%   fresthistComparison(FILES) plots 1D FRET histogram, summing over all
%   traces and frames, for each .traces file in the cell array FILES.
%
%   fresthistComparison() will prompt the user for files.
%
%   OUT = fresthistComparison(...) will return the histogram data in the matrix
%   OUT, but will not plot anything.
%
%   [...] = frethistComparison(FILES,PARAMS) specifies optional parameters:
%     'fret_axis':       bins for histogram plotting.
%     'contour_length':  number of frames to sum.
%     'pophist_offset':  number of frames to skip at the beginning.
%     'calcErrorBars':   Calculating error bars with bootstrapping.
%     'removeBlinks':    Use SKM remove dark state dwells.
%
%  frethistComparison(AX,...) plots in the the scalar axes AX.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


% Define default parameters.
% These are persistent across calls to frethistComparison, but if parameters
% are explicitly passed in (File->New/Open menus, etc), they take precedence.
persistent params;

if isempty(params)
    constants = cascadeConstants();
    params = constants.defaultMakeplotsOptions;
    params.removeBlinks  = false;
    params.calcErrorBars = false;
end

% Default colors, if there aren't too many files.
colors = [ 0      0      0    ; ...  % black
           0.75   0      0.75 ; ...  % purple
           0      0.75   0.75 ; ...  % cyan
           0      0.5    0    ; ...  % green
           0.75   0.75   0    ; ...  % yellow
           1      0      0    ; ...  % red
           0.6    0      0    ];     % dark red


%% Process input arguments
narginchk(0,3);
nargoutchk(0,3);
[varargout{1:nargout}] = deal([]);
[cax,args] = axescheck(varargin{:});

switch numel(args)
    case 0
        files = getFiles();
    case 1
        files = args{1};
    case 2
        [files,inputParams] = args{:};
        params = mergestruct(params, inputParams);
end

if ~iscell(files), files={files}; end
if numel(files)==0,  return;  end

% Only calculate histograms if hist data is requested.
if nargout>0 && isempty(cax),
    [varargout{1:nargout}] = pophist(files,params);
    return;
end

% Create a new figure if no target given.
if ~isempty(cax),
    hFig = get(cax,'Parent');
else
    hFig = figure;
end


%% Draw histograms
set(hFig,'pointer','watch'); drawnow;

nFiles = numel(files);

if nFiles>size(colors,1),
    colors = zeros(nFiles,3);
    interval = (1/(nFiles-1));
    colors(:,1) = 0:interval:1;
    colors(:,3) = 1:-interval:0;
end


% Calculate histograms
[fretaxis,frethist,errors,label] = pophist(files,params);
% output = zeros( numel(fretaxis), size(frethist,2)*2+1 );
% output(:,1) = fretaxis;
% output(:,2:2:end) = frethist;
% output(:,3:2:end) = errors;
output = [to_col(fretaxis) frethist];
% guidata(hFig, struct('files',{files},'output',output));  %for frethistComparison_save()


% Plot histograms. Splines are for better display only.
cax = newplot(hFig);
hold(cax,'on');
set(cax,'ColorOrder',colors);

sx = fretaxis(1):0.001:fretaxis(end);
sy = spline(fretaxis, frethist', sx);
plot(cax, sx, sy, 'LineWidth',3);

if params.calcErrorBars,
    for i=1:nFiles,
        errorbar(cax, fretaxis, frethist(:,i), errors(:,i)/2, '.', ...
                  'LineWidth',1, 'Color',colors(i,:));
    end
end


% Decorate the plot with axes etc.
hold(cax,'off');
ylabel( cax, 'Counts (%)' );
xlabel( cax, label );
xlim( cax, [0.1 1.0] );
yl = ylim(cax);
ylim( cax, [0 yl(2)] );

titles = trimtitles(files);
if numel(titles)>1,
    legend(cax, titles);
    title(cax,'');
else
    delete(legend(cax));
    title(cax, titles{1});
end

set(hFig,'pointer','arrow'); drawnow;



%% Add menus to change settings, get data, open new plots, etc.
fields = {'removeBlinks', 'calcErrorBars', 'contour_length', 'pophist_offset'};
prompt = {'Remove blinks:', 'Show error bars:', 'Frames', 'Offset:'};

defaultFigLayout( hFig, @(~,~)frethistComparison(getFiles(),params), ...  %File->New callback
                        @(~,~)frethistComparison(cax,getFiles(),params), ...  %File->Open callaback
                        {@exportTxt,files,output}, ...                     %Export callaback
      { 'Change settings...',@(~,~)settingdlg(params,fields,prompt,@frethistComparison,{cax,files}); ...
        'Reset settings',    @(~,~)frethistComparison(cax,files); ...
        'Copy values',      {@clipboardmat,output}  }  );

end



%%==========================================================================%%
%% Calculate 1D FRET histograms
function [fretaxis,frethist,errors,label] = pophist(files,settings)


% Display settings:
sumlen = settings.contour_length; % how many frames to use.
pophist_offset = settings.pophist_offset; % first N frames to throw out.
fretaxis = settings.fret_axis';
nbins = length(fretaxis);

if settings.calcErrorBars
    nBootstrap = 100;  %number of bootstrap samples to make.
else 
    nBootstrap = 1;
end

% Model for removing dark state noise. Adjust if any state is < 0.4.
if settings.removeBlinks
    model = QubModel(2);
    model.p0    = [0.01 0.99]';
    model.rates = [0 5; 1 0];
    model.mu       = [0.01 0.3];
    model.sigma    = [0.061 0.074];
    model.fixMu    = true(1,2);
    model.fixSigma = true(1,2);
    
    skmParams.quiet = 1;
end


%% Calculate histograms
nFiles = numel(files);
frethist = zeros(nbins,nFiles);
errors = zeros(size(frethist));
labels = cell(nFiles,1);

for i=1:nFiles
    % Load FRET data
    data = loadTraces( files{i} );
    fret = data.fret( :, pophist_offset+(1:sumlen) );
    nTraces = size(fret,1);
    labels{i} = data.fretAxisLabel;
        
    % Idealize data to 2-state model to eliminate dark-state dwells
    if settings.removeBlinks
        idl = skm( fret, data.sampling, model, skmParams );
    else
        % Use all datapoints for histogram otherwise
        idl = repmat(2,size(fret));
    end
    
    % Calculate FRET histograms from many bootstrap datasets
    pophist = zeros(nbins,nBootstrap);
    fret = fret(:,1:sumlen);
    
    for s=1:nBootstrap,
        % Construct bootstrap datasets
        if s==1,
            idxBootstrap = 1:nTraces;
        else
            idxBootstrap = floor(rand(nTraces,1)*nTraces)+1;
        end
        
        data = fret( idxBootstrap, : );     %bootstrap traces
        data = data( idl(idxBootstrap,:)==2 ); %non-zero FRET only
        
        % Create FRET histogram from the bootstrapped dataset
        histdata  = hist( data, fretaxis );
        pophist(:,s) = 100*histdata/sum(histdata);   %normalization
    end
    frethist(:,i) = pophist(:,1);
    
    % Calculate and plot error bars
    if settings.calcErrorBars
        errors(:,i) = std(pophist,[],2);
    end
end

if numel(unique(labels))==1,
    label = labels{1};
else
    label = 'FRET (various)';
end


end %function pophist




%%==========================================================================%%
%% Save histogram to file
function exportTxt(~,~,files,output)

nFiles = numel(files);

if nFiles==1,
    [p,f] = fileparts(files{1});
    outputFilename = fullfile(p, [f '_pophist.txt']);
else
    outputFilename = 'pophist.txt';
end

[f,p] = uiputfile('*.txt', [mfilename ': save histograms'], outputFilename);
outputFilename = fullfile(p,f);

if f~=0,
    % Write header line
    fid = fopen(outputFilename,'w');
    fprintf(fid,'FRET\t%s\n', strjoin(trimtitles(files),'\t'));
    fclose(fid);

    dlmwrite(outputFilename,output,'-append','delimiter','\t');
end

end



