function varargout = scaleacceptor( varargin )
% scaleacceptor Scale acceptor fluorescence so that gamma is ~1.
%
%    The FRET-distance relationship is generally defined assuming that the donor
%    and acceptor intensities have roughly the same apparent brightness
%    (gamma=1). If this is not the case, a correction is needed. This script
%    scales the acceptor fluorescence intensity by a constant factor
%
%    To scale acceptor intensity of all traces in a file:
%    
%       data = SCALEACCEPTOR( filename, scale_factor, output_filename );
%
%    If any of these parameters are not specified, you will be prompted for
%    them.
%
%    See also gammacorrect.

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
    ch = ~cellfun( @isempty, strfind(data.channelNames,'acceptor') );
    nCh = sum(ch);
    prompts = cellfun( @(x)sprintf('Scale %s by:',x), data.channelNames(ch), ...
                       'UniformOutput',false );
    
    answer = inputdlg( prompts, 'Enter factor to scale acceptor intensity', 1, ...
                                                       repmat({'1'},[nCh 1]) );
    if isempty(answer) || isempty(answer{1}),
        return;  %user hit cancel.
    end
    
    scale_factor = str2double(answer);
end



%% Scale acceptor channel and recalculate fret

% Scale each acceptor channel by the gamma factor from the user.
for i=1:numel(scale_factor),
    if i==1,
        ch = 'acceptor';
    else
        ch = sprintf('acceptor%d',i);
    end
    data.(ch) = scale_factor(i)*data.(ch);
    
    % Keep track of adjustments in trace metadata.
    idxA = strcmpi(data.channelNames,ch);
    for j=1:data.nTraces,
        data.traceMetadata(j).scaleFluor(idxA) = data.traceMetadata(j).scaleFluor(idxA).*scale_factor(i);
    end
end

% Recalculate FRET using the new acceptor fluorescence values.
data.recalculateFret();

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






