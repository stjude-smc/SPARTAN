function varargout = centerQuad( input, varargin )
%centerQuad  select traces within the center quadrant of the field of view
%
%   centerQuad(FILE) loads the .traces given file, selects traces within the
%   center quadrant of the field of view, and saves a new file with the
%   extension '_centerQuad.traces'. If there is an associated .dwt file, a new
%   '_centerQuad.qub.dwt' file will  be generated with the selected traces.
%
%   centerQuad(FILES) with a cell array of file names will process each
%   individually as described above.
%   
%   DATA = centerQuad(FILE) returns the selected traces as a Traces object and
%   will NOT save any output files.
%
%   DATA = centerQuad(DATA) modifies the Traces object DATA in place.
%
%   [DATA,IDX] = centerQuad(...) also returns the indexes of the traces that
%   were within the center quadrant.
% 
%   the movie may be excluded from analysis.

%   Copyright 2016-2017 Cornell University All Rights Reserved.


%% Process input arguments
narginchk(0,1);
nargoutchk(0,1);

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
    [varargout{1:nargout}] = centerQuad2(input, varargin{:});

elseif iscell(input)
    for i=1:numel(input)
        centerQuad2(input{i}, varargin{:});
    end %for each file
end


end %FUNCTION centerQuad



%% Function that actually remove traces at the edges
function [data,idx] = centerQuad2(input)

% Load input data
if ischar(input)
    data = loadTraces(input);
elseif isa(input,'Traces')
    data = input;
else
    error('Invalid input argument. Must be a filename or Traces object');
end

% Select only traces within the center quadrant
% ************FIXME: use mod.
x = [data.traceMetadata.donor_x];
y = [data.traceMetadata.donor_y];
stkX = max(x);
stkY = max(y);

sel = x>stkX*1/4 & x<stkX*3/4 & y>stkY*1/4 & y<stkY*3/4;
data.subset(sel);
idx = find(sel);

% Display the selected molecules to verify.
%subplot(1,numel(files),i);
%scatter( x(sel), y(sel),  'r' ); hold on;
%scatter( x(~sel),y(~sel), 'k' );
%xlim([0 stkX]); ylim([0 stkY]);

if nargout==0 && ischar(input)
    % Save the resulting traces
    [p,f,e] = fileparts(input);
    saveTraces( fullfile(p,[f '_centerQuad' e]), data );

    % If an idealization file exists, crop it also.
    dwtfname = fullfile(p,[f '.qub.dwt']);
    if ~exist(dwtfname,'file')
        dwtfname = fullfile(p,[f '.dwt']);
    end
    
    % Convert dwt to idealization, remove traces, convert to dwt, and save.
    if exist(dwtfname,'file')==2,
        [dwt,sampling,offsets,model] = loadDWT(dwtfname);
        idl = dwtToIdl(dwt,data.nFrames,offsets);
        idl = idl(sel,:);
        [dwt,offsets] = idlToDwt(idl);
        saveDWT( fullfile(p,[f '_centerQuad.qub.dwt']), dwt, offsets, model, sampling );
    end
    
end

end %FUNCTION centerQuad2



