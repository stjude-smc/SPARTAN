function varargout = crosstalkcorrect(input, varargin)
% CROSSTALKCORRECT  subtract donor to acceptor crosstalk (spectral overlap).
%
%   crosstalkcorrect( FNAME ) subtracts donor to acceptor crosstalk from all
%   traces in the .traces file given in the string FNAME, using a single fixed
%   value that is the median of apparent crosstalk values in the traces.
%   The corrected traces are saved to a new file with the extension '_ccor'.
%
%   crosstalkcorrect( FILES ) corrects each file in the cell array FILES,
%   calculating a distinct crosstalk value for each file. 
%
%   [DATA,CT] = crosstalkcorrect( FNAME ) returns the corrected data and
%   crosstalk value used for correction instead of saving to file.
%
%   DATA = crosstalkcorrect( DATA ) corrects the TracesFret object DATA.
%   (NOTE: input is modified in place. Create a copy to preserve the original.)
%
%   [...] = crosstalkcorrect( ..., CT ) uses a user-supplied crosstalk value
%   instead of calculating one.
%
%   ALGORITHM:
%   Spectral crosstalk, where a fraction of donor fluorescence appears on the
%   acceptor channel, is estimated from the elevated acceptor channel baseline
%   after acceptor photobleaching. Bleaching is determined as the last point
%   above 0.2 FRET. This algorithm will not work properly if crosstalk is
%   > 15% or if there are non-zero FRET states with low efficiency also in that
%   range. In such cases, manually correct with an estimated value first.
%
%   See also: calc_crosstalk, scaleacceptor, gammacorrect, correctTraces.

%   Copyright 2014-2017 Cornell University All Rights Reserved.


% Check input arguments
narginchk(0,2);
nargoutchk(0,2);

if nargin<1 || isempty(input),
    % If no files given, obtain a list from the user.
    input = getFiles('*.traces;*.rawtraces');
    if isempty(input), return; end  %user hit cancel
    
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
    [varargout{1:nargout}] = ccorr(input, varargin{:});

elseif iscell(input)
    for i=1:numel(input)
        ccorr(input{i}, varargin{:});
    end %for each file
end


end %function crosstalkcorrect



%% Actually correct the crosstalk
function varargout = ccorr(input, ctval)

if ischar(input)
    data = loadTraces(input);
elseif isa(input,'Traces')
    data = input;
else
    error('Invalid input argument. Must be a filename or Traces object');
end

if isChannel(data,'acceptor2'),
    disp('Warning: this function may not work correctly with multi-color FRET');
end

% Estimate the crosstalk if not given by the user.
if nargin<2
    ctval = calc_crosstalk(data);
    disp(ctval);
end

% Add values to existing crosstalk matrix
idxD  = strcmpi(data.channelNames,'donor');
idxA1 = strcmpi(data.channelNames,'acceptor');
crosstalk = cat(3, data.traceMetadata.crosstalk);
crosstalk(idxD,idxA1,:) = crosstalk(idxD,idxA1,:) + ctval;

% Make the crosstalk correction and recalculate FRET
data = correctTraces(data, crosstalk);
data.recalculateFret();

% Save resulting data
if ischar(input) && nargout==0
    [p,f,e] = fileparts( input );
    outFilename = fullfile(p, [f '_ccorr' e]);

%     if nFiles==1,
%         [f,p] = uiputfile(outFilename,'Save corrected file');
%         if ~ischar(f), return; end
%         outFilename = fullfile(p,f);
%     end

    saveTraces( outFilename, data );
end

output = {data,ctval};
[varargout{1:nargout}] = output{1:nargout};

end %FUNCTION CCORR


