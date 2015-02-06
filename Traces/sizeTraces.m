function [nTraces,nFrames,channelNames] = sizeTraces( filenames, dim )
% Given a file name or cell array of filenames, load only the header of traces
% files to determine the 
%
% See also loadTraces.m and saveTraces.m
%

% If not given, ask the user for filenames
if nargin<1 || isempty(filenames),
    filenames = getFiles( {'*.traces;*.rawtraces','Binary Traces Files (*.traces)'; ...
              '*.traces','Traces Files (*.traces)'; ...
              '*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
              '*.*','All Files (*.*)'} );
end
if ischar(filenames),
    filenames = {filenames};
end

nFiles = numel(filenames);
if nFiles==0,
    return;
end

nTraces = zeros(nFiles,1);
nFrames = zeros(nFiles,1);
channelNames = cell(nFiles,1);


%%
for i=1:numel(filenames)
    % 1) Open the traces file.
    fid = fopen(filenames{i},'r');

    % 2) Read header information
    if fread(fid,1,'*uint32')~=0,
        error('Traces format version not supported or invalid file');
    end

    % Check validity of header data.
    magic     = fread( fid, [1,4], '*char' );  %format identifier ("magic")
    version   = fread( fid, 1, '*uint16' );    %format version number

    assert( strcmp(magic,'TRCS'), 'Invalid traces file' );
    assert( version>=3 && version<=4, 'Traces format version not supported!' );

    fread( fid, 1, '*uint8' );  %data type/precision
    fread( fid, 1, '*uint8' );  %number of channels
    nTraces(i) = fread( fid, 1, 'uint32' );
    nFrames(i) = fread( fid, 1, 'uint32' );

    % 3) Read data channel names (version 4+)
    if version>3 && nargout>=3,
        szNames = fread( fid, 1, 'uint32' );
        ch = strtrim(  fread( fid, [1,szNames], '*char' )  );
        c = textscan(ch, '%s', 'Delimiter',char(31));
        ch = c{1}';

        % Remove empty (trailing) channel names. This might happen if extra
        % delimiters are added to the end to pad to word boundries.
        channelNames{i} = ch( ~cellfun(@isempty,ch) );
    end
end

% Allow the user to specify the dimension to get, like size().
if nargin>1,
    assert( isnumeric(dim) && numel(dim)==1 && dim>=1 && dim<=2, 'Invalid size dimension' );
    if dim==2,
        nTraces=nFrames;
    end
end


end %function tracesFileSize
