function [viewer,movieFilename] = showMovie(varargin)
%showMovie  Display wide-field movie that produced a particular trace.
%
%   AX = showMovie(DATA, MOL) displays the raw data movie associated with a
%   molecule number (MOL) from the Traces object DATA in a new window
%   (The movie filename is determined from traceMetadata.ids).
%   The selected molecule will also be highlighted in all spectral channels.
%
%   showMovie(AX, DATA, MOL) updates an existing display to highlight a
%   different molecule number (MOL).
%
%   PARAMS struct must contain: geometry, wavelengths, chNames, chDesc, 
%     nAvgFrames, bgBlurSize (see @MovieParser/openStk.m for details).
%
%   See also: sorttraces, MovieViewer.

%   Copyright 2017-2018 Cornell University All Rights Reserved.


% Process input arguments.
narginchk(3,4); nargoutchk(0,2);
viewer = [];

if nargin>=1 && isa(varargin{1},'MovieViewer')
    viewer = varargin{1};
    [data,m,fpath] = varargin{2:end};
else
    [data,m,fpath] = varargin{:};
end


% Parse trace identifiers to predict original movie filename.
id = data.traceMetadata(m).ids;
assert( any(id=='#') );
output = strsplit(id,'#');
[moviePath,~] = deal( output{:} );

% Look for the movie in the directory containing the loaded traces file.
% using .tif extension if it is .rawtraces etc (in some old versions).
[~,f,e] = fileparts2(moviePath);
if ~ismember(e, {'.stk','.tif','.tiff','.pma'}), e='.tif'; end
movieFilename = fullfile(fpath, [f e]);

if ~exist( movieFilename, 'file' )
    error('Corresponding movie file not found');
end

% If trace is from a different file than the one currently loaded,
% close the existing viewer; user will have to open the new one.
% FIXME: would be nice if the transition were seamless, but this is easier.
if ~isempty(viewer) && isvalid(viewer)
    [~,fold] = fileparts2(viewer.chExtractor.movie.filename{1});
    if ~strcmp(f, fold)
        delete(viewer);
        return;
    end
end

% If the viewer is not open, launch a new one.
if isempty(viewer) || ~isvalid(viewer)
    
    % Assemble movie parsing parameters.
    try
        % Use trace metadata to create a ChannelExtractor
        ch = struct( 'name',data.fileMetadata.chDesc, ...
                     'wavelength',num2cell(data.fileMetadata.wavelengths) );
        ex = ChannelExtractor( movieFilename, data.fileMetadata.geometry, ch );
        viewer = MovieViewer( ex );
    catch e
        % If geometry field is not available (older file), use movie file
        % metadata or default to single-color if not present.
        disp('Using trace metadata to create ChannelExtractor failed; using default settings instead');
        disp(e.message);
        viewer = MovieViewer( movieFilename );
    end

    % Show movie viewer window
    viewer.show();
end


% Draw circle around the currently-selected molecule.
fluorCh = data.channelNames(data.idxFluor);
coord = cell( numel(fluorCh), 1 );

for i=1:numel(fluorCh)
    try
        x = data.traceMetadata(m).( [fluorCh{i} '_x'] );
        y = data.traceMetadata(m).( [fluorCh{i} '_y'] );
        coord{i} = [x y];
    catch e
        disp(e.message);
    end
end
viewer.highlightMolecule(coord);


end %function showMovie



function [p,f,e] = fileparts2(fname)
% Extract path, file name, and extension of a file.
% The built-in fileparts only separates paths using the current system's
% filesep, which may be different from the one used when the file was saved.

pidx = find( fname=='\'|fname=='/', 1, 'last' );
eidx = find( fname=='.', 1, 'last' );

p = fname(1:pidx);
f = fname(pidx+1:eidx-1);
e = fname(eidx:end);

end


