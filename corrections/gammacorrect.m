function varargout = gammacorrect(input, varargin)
% gammacorrect scale acceptor intensity by estimated gamma values.
%
%   gammacorrect( FNAME ) Corrects for unequal sensitivity and/or quantum yield
%   ('apparent brightness') of donor and acceptor fluorophores in two-color
%   FRET traces by scaling the donor intensity. This scaling factor (gamma) is
%   estimated for each trace as the magnitude increase in donor intensity
%   divided by the magnitude decrease in acceptor intensity upon acceptor
%   photobleaching. A single fixed gamma value is used to make the correction,
%   obtained from the median of apparent trace gamma values, excluding outliers.
%   FNAME is the path of the .traces file to correct. The corrected traces are
%   automatically saved in a new file with the extension '_gcorr.traces'.
%
%     NOTE: for best results, use this function after crosstalk and other
%     corrections and trace selection in autotrace.
%   
%     See Ha et al (1999) PNAS 95, p.893 for original algorithm description.
%
%   gammacorrect( FILES ) corrects each file in the cell array FILES,
%   calculating a distinct gamma value for each file.
%
%   [DATA,GAMMA] = gammacorrect( FNAME ) returns the corrected data and
%   gamma value used for correction instead of saving to file.
%
%   DATA = gammacorrect( DATA ) corrects the TracesFret object DATA directly.
%   (NOTE: input is modified in place. To preserve the original, create a copy.)
%
%   [...] = crosstalkcorrect( ..., GAMMA ) uses a user-supplied gamma value
%   instead of calculating one.
% 
%   See also: calc_gamma, scaleacceptor, crosstalkcorrect, correctTraces.

%   Copyright 2014-2017 Cornell University All Rights Reserved.


%% Process input arguments
narginchk(0,2);
nargoutchk(0,2);

if nargin<1 || isempty(input),
    % If no files given, obtain a list from the user.
    filter = {'*.*traces','Traces files (*.traces,*.rawtraces)';
              '*.*','All files (*.*)'};
    input = getFiles(filter);
    if isempty(input), return; end  %user hit cancel.
    
elseif iscell(input)
    if ~all(cellfun(@ischar,input))
        error('Input must be a Traces object, filename string, cell array of filenames');
    end
    if nargout>0
        error('Output only supported for scalar input')
    end
end

% Correct traces and save/return result
if isa(input,'Traces') || ischar(input)
    [varargout{1:nargout}] = gcorr(input, varargin{:});

elseif iscell(input)
    for i=1:numel(input)
        gcorr(input{i}, varargin{:});
    end
end


end %function gammacorrect





%% Actually correct for gamma.
function varargout = gcorr(input, gamma)

% Load fluorescence data
if ischar(input)
    data = loadTraces( input );
elseif isa(input,'TracesFret')
    data = input;
else
    error('Invalid input. Must be TracesFret object or path to .traces file');
end

if isChannel(data,'acceptor2'),
    disp('Warning: this function may not work correctly with multi-color FRET');
end

% Estimate the mean gamma value across all traces in the file
if nargin<2
    gamma = calc_gamma(data);
end

% Keep track of adjustments in trace metadata.
idxA1 = strcmpi(data.channelNames,'acceptor');
[data.traceMetadata.scaleFluor] = to_col(data.traceMetadata.scaleFluor);
scaleFluor = cat(2,data.traceMetadata.scaleFluor);
scaleFluor(idxA1,:) = scaleFluor(idxA1,:)*gamma;

data = correctTraces(data,[],scaleFluor);
data.recalculateFret();

% Save resulting data
if ischar(input) && nargout==0
    [p,f,e] = fileparts( input );
    outFilename = fullfile(p, [f '_gcorr' e]);

    %     if nFiles==1,
    %         [f,p] = uiputfile(outFilename,'Save corrected file');
    %         if ~ischar(f), return; end
    %         outFilename = fullfile(p,f);
    %     end

    saveTraces( outFilename, data );
end

output = {data,gamma};
[varargout{1:nargout}] = output{1:nargout};

end




