function [nTraces,nFrames,channelNames] = sizeTraces( filenames, dim )
%sizeTraces   Number and length of traces in a file
%
%   [N,M,C] = sizeTraces( FILENAME )
%   returns the number of traces (N), trace length in frames (M), and names of
%   data channels (C, a cell array of strings) from the specified FILENAME
%   (a .traces file). Only the header is loaded, which is very fast.
%
%   [N,M,C] = sizeTraces( FILENAMES )
%   returns a column vector of values, one for each cell in FILENAMES.
%
%   X = sizeTraces( FILENAME, DIM )
%   returns only the desired dimension (1=N, 2=M, 3=C).
%
%   See also loadTraces.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% If not given, ask the user for filenames
if nargin<1 || isempty(filenames),
    filenames = getFiles( {'*.traces;*.rawtraces','All Traces Files (*.traces)'; ...
              '*.traces','Traces Files (*.traces)'; ...
              '*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
              '*.*','All Files (*.*)'} );
end
if ischar(filenames),
    filenames = {filenames};
end

nFiles = numel(filenames);
nTraces = zeros(nFiles,1);
nFrames = zeros(nFiles,1);
channelNames = cell(nFiles,1);


%%
for i=1:numel(filenames)
    % If loading a .txt file, we have no choice but to load the full file.
    [~,~,e] = fileparts(filenames{i});
    if strcmp(e,'.txt'),
        data = loadTraces(filenames{i});
        nTraces(i) = data.nTraces;
        nFrames(i) = data.nFrames;
        channelNames{i} = {'donor','acceptor','fret'};
        continue;
    end
    
    % 1) Open the traces file.
    fid = fopen(filenames{i},'r');

    % Handle legacy (pre 2012) format .traces files
    if fread(fid,1,'*uint32')~=0,
        frewind(fid);
        nFrames(i) = fread(fid,1,'int32');
        nTraces(i) = fread(fid,1,'int16')/2;  %D,A are counted separately
        channelNames{i} = {'donor','acceptor'};
        fclose(fid);
        continue;
    end

    % 2) Read header information.
    % Check validity of header data.
    magic     = fread( fid, [1,4], '*char' );  %format identifier ("magic")
    version   = fread( fid, 1, '*uint16' );    %format version number

    assert( strcmp(magic,'TRCS'), 'Invalid traces file' );
    assert( version>=3 && version<=5, 'Traces format version not supported!' );

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
    
    fclose(fid);
end

% Allow the user to specify the dimension to get, like size().
if nargin>1,
    assert( isnumeric(dim) && numel(dim)==1 && dim>=1 && dim<=2, 'Invalid size dimension' );
    if dim==2,
        nTraces=nFrames;
    end
end


end %function tracesFileSize
