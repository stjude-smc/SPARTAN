function varargout = frethistComparison(files, inputParams)
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

%   Copyright 2007-2016 Cornell University All Rights Reserved.


% Define default parameters.
constants = cascadeConstants();
params = constants.defaultMakeplotsOptions;
params.removeBlinks  = false;
params.calcErrorBars = false;

if nargin>1,
    params = mergeStruct(params, inputParams);
end

% Prompt user for filenames if not supplied
if nargin<1,
    files = getFiles();
end
if ~iscell(files), files = {files}; end
if numel(files)==0,  return;  end


% Only calculate histograms if hist data is requested.
if nargout>0,
    [fretaxis,frethist,errors] = frethistComparison_calc(files,params);
    output = {fretaxis,frethist,errors};
    [varargout{1:nargout}] = output{1:nargout};
    return;
end


% Save the histogram data in the figure for access by callback functions.
handles = struct('hFig',figure, 'files',{files}, 'params',params);
guidata(handles.hFig,handles);

% Add menu items
hTxtMenu = findall(gcf, 'tag', 'figMenuGenerateCode');
set(hTxtMenu, 'Label','Export as .txt', 'Callback',@frethistComparison_save);

hEditMenu = findall(gcf, 'tag', 'figMenuEdit');
delete(allchild(hEditMenu));
uimenu('Label','Display settings...', 'Parent',hEditMenu, 'Callback',@frethistComparison_settings);

% Plot the histograms
frethistComparison_display(handles.hFig);


end %function frethistComparison



%%
function [fretaxis,frethist,errors] = frethistComparison_calc(files,settings)
%


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
if settings.removeBlinks,
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
nFiles = numel(files);
frethist = zeros(nbins,nFiles);
errors = zeros(size(frethist));

for i=1:nFiles
    % Load FRET data
    data = loadTraces( files{i} );
    fret = data.fret( :, pophist_offset+(1:sumlen) );
    [nTraces,nFrames] = size(fret);
        
    % Idealize data to 2-state model to eliminate dark-state dwells
    if settings.removeBlinks,
        [dwt,~,~,offsets] = skm( fret, data.sampling, model, skmParams );
        idl = dwtToIdl( dwt, offsets, nFrames, nTraces );
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


end %function frethistComparison_calc





%% Display histograms
function frethistComparison_display(hObject,~,~)
%

handles = guidata(hObject);
set(handles.hFig,'pointer','watch'); drawnow;

params  = handles.params;
titles  = trimtitles(handles.files);
nFiles  = numel(titles);


% Define parameters
colors = [ 0      0      0    ; ...  % black
           0.75   0      0.75 ; ...  % purple
           0      0.75   0.75 ; ...  % cyan
           0      0.5    0    ; ...  % green
           0.75   0.75   0    ; ...  % yellow
           1      0      0    ; ...  % red
           0.6    0      0    ];     % dark red

if nFiles>size(colors,1),
    colors = zeros(nFiles,3);
    interval = (1/(nFiles-1));
    colors(:,1) = 0:interval:1;
    colors(:,3) = 1:-interval:0;
end


% Calculate histograms anew
[fretaxis,frethist,errors] = frethistComparison_calc(handles.files,params);
% handles.output = zeros( numel(fretaxis), size(frethist,2)*2+1 );
% handles.output(:,1) = fretaxis;
% handles.output(:,2:2:end) = frethist;
% handles.output(:,3:2:end) = errors;
handles.output = [to_col(fretaxis) frethist];
guidata(hObject,handles);


% Plot histograms. Splines are for better display only.
cax = cla(handles.hFig);
hold(cax,'on');
set(cax,'ColorOrder',colors);

sx = fretaxis(1):0.001:fretaxis(end);
sy = spline( fretaxis, frethist', sx );
plot( cax, sx, sy, 'LineWidth',3 );

if params.calcErrorBars,
    for i=1:nFiles,
        errorbar( fretaxis, frethist(:,i), errors(:,i)/2, '.', ...
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

set(handles.hFig,'pointer','arrow'); drawnow;


end





%% Dialog to change parameters
function frethistComparison_settings(hObject,~,~)

handles = guidata(hObject);
opt = handles.params;

% 1. Get the new value from the user.
prompt = {'Remove blinks:', 'Show error bars:', 'Frames', 'Offset:'};
fields = {'removeBlinks', 'calcErrorBars', 'contour_length', 'pophist_offset'};
currentopt = cellfun( @(x)num2str(opt.(x)), fields, 'UniformOutput',false );

answer = inputdlg(prompt, [mfilename ' display settings'], 1, currentopt);
if isempty(answer), return; end  %user hit cancel

% 2. Save new parameter values from user.
for k=1:numel(answer),
    original = opt.(fields{k});
    
    if isnumeric(original)
        opt.(fields{k}) = str2double(answer{k});
    elseif islogical(original)
        opt.(fields{k}) = logical(str2double(answer{k}));
    else
        opt.(fields{k}) = answer{k};
    end
    
    % FIXME: verify string fields.
end

handles.params = opt;
guidata(hObject,handles);

% 3. Redraw plots.
frethistComparison_display(hObject);



end %function dwellhist_settings





%% Save histogram to file
function frethistComparison_save(hObject,~,~)

handles = guidata(hObject);
files = handles.files;
output = handles.output;
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



