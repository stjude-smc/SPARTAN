function varargout = avgFretTime(varargin)
%avgFretTime  Average FRET trajectory across all traces with bleaching removed.
%
%   OUT = avgFretTime(FILES) median FRET trajectory for each .traces file in
%   the cell array FILES (first column is time axis in seconds).
%
%   OUT = avgFretTime() will prompt the user for files to load.
%
%   avgFretTime(...) with no outputs displays the data in a new figure.
%
%   avgFretTime(AX,...) plots in the the scalar axes AX.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


% Default parameter values
params.truncateLen = 300;  %frames to calculate over
params.min_fret = 0.175;  % minimum fret value, below which we assume there is no FRET.



%% Process input arguments
narginchk(0,3);
nargoutchk(0,1);
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
nFiles = numel(files);
if nFiles==0,  return;  end



%% Calculate average FRET trajectories
output = zeros(params.truncateLen,1+nFiles);

for i=1:numel(files),
    % Load FRET data and truncate to target length
    data = loadTraces( files{i} );
    data.nFrames = params.truncateLen;
    
    if i==1,
        output(:,1) = data.time;  %convert to seconds
    else
        if ~all(data.time==output(:,1)'),
            warning('Time axes do not match!');
        end
    end
    
    % For each trace, average the FRET values at each frame to create an
    % average FRET trace. Exclude dark states and extreme values.
    for j=1:params.truncateLen,
        nonzero = data.fret(:,j) >= params.min_fret & data.fret(:,j)<1.2;
        output(j,i+1) = mean( data.fret(nonzero,j) );
    end

end %for each trace.

output(:,1) = output(:,1)/1000;  %convert to seconds


if nargout>0 && isempty(cax),
    varargout{1} = output;
    return;
end


%% Plot the result
if ~isempty(cax),
    hFig = get(cax,'Parent');
else
    hFig = figure;
end
cax = newplot(hFig);

plot( cax, output(:,1), output(:,2:end) );
xlabel(cax, 'Time (s)');
ylabel(cax, 'Average FRET value');

titles = trimtitles(files);
if nFiles>1,
    title('');
    legend(cax, titles);
else
    title(cax, titles{1});
    delete(legend(cax));
end


%% Add menu items
defaultFigLayout( hFig, @(~,~)avgFretTime(getFiles(),params), ...
                        @(~,~)avgFretTime(cax,getFiles(),params),  ...
      {@exportTxt,files,output}, {'Copy values',{@clipboardmat,output}} );


end %function avgFretTime




%  =========================================================================  %
%% Save the output to file.
function exportTxt(~,~,files,output)

if numel(files)==1,
    [p,f] = fileparts(files{1});
    outputFilename = fullfile(p, [f '_' mfilename '.txt']);
else
    outputFilename = [mfilename '.txt'];
end

[f,p] = uiputfile('*.txt', [mfilename ': save output'], outputFilename);
outputFilename = fullfile(p,f);

if f~=0,
    % Write header line
    fid = fopen(outputFilename,'w');
    fprintf(fid,'Time (s)\t%s\n', strjoin(trimtitles(files),'\t'));
    fclose(fid);

    dlmwrite(outputFilename,output,'-append','delimiter','\t');
end

end %function avgFretTime_save



