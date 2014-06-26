function data = scale_acceptor( varargin )
% SCALE_ACCEPTOR:  Scale acceptor fluorescence so that gamma is ~1.
%
%    The FRET-distance relationship is generally defined assuming that the donor
%    and acceptor intensities have roughly the same apparent brightness
%    (gamma=1). If this is not the case, a correction is needed. This script
%    scales the acceptor fluorescence intensity by a constant factor
%
%    To scale acceptor intensity of all traces in a file:
%    
%       data = SCALE_ACCEPTOR( filename, scale_factor, output_filename );
%
%    If any of these parameters are not specified, you will be prompted for
%    them.
%
%    See also GAMMA_CORRECT.
%


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

% Get scale factor. If none given, use a default of 5, which is common for
% Cy3/Cy7 FRET, which is the main reason you would use this script.
if nargin>=2,
    scale_factor = varargin{2};
else
    answer = inputdlg( 'Enter factor by which to scale the acceptor:', ...
                       'scale_acceptor: enter scale factor', 1, {'5'} );
    if isempty(answer) || isempty(answer{1}),
        return;  %user hit cancel.
    end
    
    scale_factor = str2double(answer);
end



%% Scale acceptor channel and recalculate fret
data.acceptor = scale_factor*data.acceptor;

fret = data.acceptor ./ (data.acceptor+data.donor);
fret( data.fret==0 ) = 0;  %copy marker for when donor is dark from original data.
data.fret = fret;

clear fret;


%% Save the result

% If no output filename given, generate one.
if nargin>=3,
    out_filename = varargin{3};
else
    [p,f] = fileparts(filename);    
    [f,p] = uiputfile( '*.traces', 'Select output filename', fullfile(p,f) );
    if ischar(f),
        out_filename = fullfile(p,f);
    else
        return;  %user hit cancel.
    end
end

% Save the output.
saveTraces( out_filename, data );





end %function scale_acceptor






