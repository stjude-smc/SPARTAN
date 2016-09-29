function varargout = scaleacceptor( varargin )
% scaleacceptor Scale acceptor fluorescence so that gamma is ~1.
%
%    Multiplies fluorescence channels by a user-specified factor to correct
%    for unequal brightness/detection efficiency. This correction is
%    cumulative to any previous corrections made:
%    
%       data = SCALEACCEPTOR( filename, scale_factor, output_filename );
%
%    If any of these parameters are not specified, you will be prompted for
%    them.
%
%    See also gammacorrect, correctTraces.

%   Copyright 2014-2016 Cornell University All Rights Reserved.


%% Process input

% Get filename and load fluorescence data.
if nargin<1 || isempty(varargin{1}),
    % If no file is specified, ask for one from the user.
    [f,p] = uigetfile( {'*.traces;*.rawtraces','Binary Traces Files (*.traces;*.rawtraces)'; ...
                        '*.txt','Old format traces files (*.txt)'; ...
                        '*.*','All Files (*.*)'}, 'Select a traces file');
    if p==0, return; end  %user hit cancel.
    filename = fullfile(p,f);
    data = loadTraces(filename);

elseif nargin>=1 && ischar(varargin{1}),
    filename = varargin{1};
    data = loadTraces(filename);

else
    % Otherwise, assume we are given a data struct or Traces object.
    assert( isstruct(varargin{1}) || isa(varargin{1},'Traces') );
    filename = '';
    data = varargin{1};
end

% Get scale factor from the user if not given on the command prompt.
if nargin>=2,
    scale_factor = varargin{2};
else
    % Get current value for reference
    scale_factor = mean( cat(2, to_col(data.traceMetadata.scaleFluor)), 2 );
    defaults = cellfun( @num2str, num2cell(scale_factor), 'Uniform',false );
    
    prompts = cellfun( @(x)sprintf('Scale %s by:',x), ...
                   data.channelNames(data.idxFluor), 'Uniform',false );
    
    answer = inputdlg( prompts, 'Enter factor to scale intensities', ...
                                                               1, defaults);
    if isempty(answer) || isempty(answer{1}),
        return;  %user hit cancel.
    end
    
    scale_factor = str2double(answer);
end



%% Scale acceptor channel and recalculate fret
data = correctTraces(data, [], scale_factor);
data.recalculateFret();

% FIXME: if return requested and no filename, should just return data.
if nargout>=1,
    varargout{1} = data;
end


%% Save the result

% If no output filename given, generate one.
if nargin>=3,
    out_filename = varargin{3};
else
    [p,f,e] = fileparts(filename);
    [f,p] = uiputfile( '*.traces', 'Select output filename', fullfile(p,[f '_corr' e]) );
    if ischar(f),
        out_filename = fullfile(p,f);
    else
        return;  %user hit cancel.
    end
end

% Save the output.
saveTraces( out_filename, data );





end %function scaleacceptor






