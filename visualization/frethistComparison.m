function varargout = frethistComparison(files, settings)
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
%   [...] = frethistComparison(FILES,PARAMS) specifies optional parameters in
%   the struct PARAMS. Many are the same as for makeplots, but also:
%
%     'removeDarkState': Use SKM to idealize FRET to detect dark state dwells
%                        (acceptor blinking) and remove from histograms.
%
%     'model':           Model for idealizing FRET data to remove dark states.
%
%     'nBootstrap':      Number of bootstrap samples for calculating error bars.
%                        If set to 1, no error bars will be shown.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


if nargin<3,
    constants = cascadeConstants();
    settings = constants.defaultMakeplotsOptions;
end

% Display settings:
pophist_sumlen = settings.contour_length; % how many frames to use.
pophist_offset = settings.pophist_offset; % first N frames to throw out.
fretaxis = settings.fret_axis';
nbins = length(fretaxis);

% Settings for removing the zero-FRET state.
settings.removeDarkState = true;

% Settings for error bar calculation.
settings.calcErrorBars = true;

if settings.calcErrorBars
    nBootstrap = 100;  %number of bootstrap samples to make.
else 
    nBootstrap = 1;
end


% Prompt user for filenames if not supplied
if nargin<1,
    files = getFiles();
end
if ~iscell(files), files = {files}; end

nFiles = numel(files);
if nFiles==0,  return;  end


% Model for removing dark state noise. Adjust if any state is < 0.4.
if settings.removeDarkState,
    model.p0    = [0.01 0.99]';
    model.rates = [0 5; 1 0];
    model.mu       = [0.01 0.3];
    model.sigma    = [0.061 0.074];
    model.fixMu    = true(1,2);
    model.fixSigma = true(1,2);
    model = QubModel(model);
    
    skmParams.quiet = 1;
end


%% Calculate histograms
frethist = zeros(nbins,2*nFiles);  %hist1, err1, hist2, err2, ...

for i=1:nFiles
    
    % Load FRET data
    data = loadTraces( files{i} );
    fret = data.fret( :, pophist_offset+(1:pophist_sumlen) );
    [nTraces,nFrames] = size(fret);
        
    % Idealize data to 2-state model to eliminate dark-state dwells
    if settings.removeDarkState,
        [dwt,~,~,offsets] = skm( fret, data.sampling, model, skmParams );
        idl = dwtToIdl( dwt, offsets, nFrames, nTraces );
    else
        % Use all datapoints for histogram otherwise
        idl = repmat(2,size(fret));
    end
    
    % Calculate FRET histograms from many bootstrap datasets
    pophist = zeros(nbins,nBootstrap);
    fret = fret(:,1:pophist_sumlen);
    
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
    
    % Calculate and plot error bars
    if settings.calcErrorBars
        pophistErrors = std(pophist,[],2);
        frethist(:,2*i) = pophistErrors;
    end
    
    % Add histogram from current dataset to output
    % Add bootstrapped errors from current dataset
    frethist(:,2*i-1) = pophist(:,1); %use the first set: all traces.
end


% Save the results. If calcErrorBars is true, every other column has the
% error bars of the associated histogram.
if ~settings.calcErrorBars,
    frethist = frethist(:,1:2:end); %remove error bar columns.
end


output = [fretaxis frethist];
if nargout>0,
    varargout{1} = output;
    return;
end

%% Setup GUI

% Create titles if not specified.
if nargin<2,
    titles = trimtitles(files);
end

hFig = figure;

% Save the histogram data in the figure for access by callback functions.
handles.hFig = hFig;
handles.titles = titles;
handles.files = files;
handles.params = settings;
handles.frethist = frethist;
handles.fretaxis = fretaxis;
handles.output = output;
guidata(hFig,handles);

% Add menu items
hTxtMenu = findall(gcf, 'tag', 'figMenuGenerateCode');
set(hTxtMenu, 'Label','Export as .txt', 'Callback',@frethistComparison_save);

hEditMenu = findall(gcf, 'tag', 'figMenuEdit');
delete(allchild(hEditMenu));
uimenu('Label','Display settings...', 'Parent',hEditMenu, 'Callback',@frethistComparison_display);


frethistComparison_display(hFig,handles);


end %function frethistComparison




%% Display histograms
function frethistComparison_display(hObject,~,~)


colors = [ 0      0      0    ; ...  % black
           0.75   0      0.75 ; ...  % purple
           0      0.75   0.75 ; ...  % cyan
           0      0.5    0    ; ...  % green
           0.75   0.75   0    ; ...  % yellow
           1      0      0    ; ...  % red
           0.6    0      0    ];     % dark red


handles = guidata(hObject);
frethist = handles.frethist;
fretaxis = handles.fretaxis;
titles   = handles.titles;
params   = handles.params;
nFiles = numel(titles);

cax = cla(handles.hFig);
hold(cax,'on');

% If there are a ton of datapoints, we can't use the simple colors above.
% Instead, just use a simple blue-to-red gradient
if nFiles>size(colors,1),
    colors = zeros(nFiles,3);
    interval = (1/(nFiles-1));
    colors(:,1) = 0:interval:1;
    colors(:,3) = 1:-interval:0;
end
set(cax,'ColorOrder',colors);

% Spline interpolate the data so it's easier to see (but not saved that way)
sx = fretaxis(1):0.001:fretaxis(end);
sy = spline( fretaxis, frethist(:,2*(1:nFiles)-1)', sx );
plot( cax, sx, sy, 'LineWidth',3 );

if params.calcErrorBars,
    for i=1:nFiles,
        errorbar( fretaxis, frethist(:,2*i-1), frethist(:,2*i)/2, '.', ...
                  'LineWidth',1, 'Color',colors(i,:) );
    end
end


% Decorate the plot with axes etc.
hold(cax,'off');
ylabel( cax, 'Counts (%)' );
xlabel( cax, 'FRET' );
xlim( cax, [0.1 1.0] );
yl = ylim(cax);
ylim( cax, [0 yl(2)] );

if nFiles>1,
    legend(cax, titles);
else
    title(cax, titles{1});
end



end



%% Save histogram to file
function frethistComparison_save(hObject,~,~)

handles = guidata(hObject);
files = handles.files;
output = handles.output;
titles = handles.titles;
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
    fprintf(fid,'FRET\t%s\n', strjoin(titles,'\t'));
    fclose(fid);

    dlmwrite(outputFilename,output,'-append','delimiter','\t');
end

end



