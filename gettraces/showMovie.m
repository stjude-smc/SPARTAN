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
if ~ismember(e, {'.stk','.tif','.tiff'}), e='.tif'; end
movieFilename = fullfile(fpath, [f e]);

if ~exist( movieFilename, 'file' )
    error('Corresponding movie file not found');
end


% The viewer has been closed or never opened, create one.
if isempty(viewer) || ~isvalid(viewer)
    
    % Assemble movie parsing parameters.
    constants = cascadeConstants;
    params = constants.gettraces_profiles(1);  %single color
    params.chNames = data.channelNames;
    try
        params.wavelengths = data.fileMetadata.wavelengths;
        params.chDesc      = data.fileMetadata.chDesc;
        params.geometry    = data.fileMetadata.geometry;  %new in v3.7
    catch e
        % If geometry field is not available (older file), gracefully fail by
        % defaulting to single-color mode (one axis instead of three).
        disp(e.message);
    end

    % Create and show movie viewer
    viewer = MovieViewer( movieFilename, params );
    viewer.show();
    
    
% If trace is from a different file than the one currently load,
% load the new movie quietly into the existing window.
else
    [~,fold] = fileparts2(viewer.movie.filename{1});
    if ~strcmp(f, fold),
        viewer.load(movieFilename);
    end
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
viewer.highlightPeaks(coord);


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


