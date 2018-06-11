function [axFOV,movieFilename] = showMovie(varargin)
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
%   See also: sorttraces, gettraces.

%   Copyright 2017 Cornell University All Rights Reserved.

% TODO: split fields; basically a small subset of gettraces functionality.
% This requires knowing the geometry, which is not available in metadata...


narginchk(2,3); nargoutchk(0,2);
[axFOV,args] = axescheck(varargin{:});
[data,m] = args{:};


%% Step 1: parse trace identifiers to predict original movie filename.
try
    id = data.traceMetadata(m).ids;
    assert( any(id=='#') );
    output = strsplit(id,'#');
    [traceID,~] = deal( output{:} );
catch
    disp('Failed to find movie source file. Metadata not available.');
    return;
end


%% If not already drawn, show movie image in new window.
if isempty(axFOV) || ~ishandle(axFOV)
    
    % Look for the movie in the current directory,
    % using .tif extension if it is .rawtraces etc (in some old versions).
    [~,f,e] = fileparts2(traceID);
    if ~ismember(e, {'.stk','.tif','.tiff'}), e='.tif'; end
    movieFilename = fullfile(pwd, [f e]);
    
    % If not found in the current folder, ask the user.
    if ~exist( movieFilename, 'file' )
        [f,p] = uigetfile('*.tif;*.tiff;*.stk', 'Select movie', movieFilename);
        if isequal(f,0), return; end  %user hit cancel
        movieFilename = fullfile(p,f);
    end

    % Create a viewer to display movie
    viewer = MovieViewer( {movieFilename} );
    axFOV = viewer.show();
end

setappdata(axFOV,'traceID',traceID);


%% Step 3. Draw circle around the currently-selected molecule.
fluorCh = data.channelNames(data.idxFluor);
nCh = numel(fluorCh);

% Get x/y location of current molecule in each fluorescence channel.
coord = nan(nCh,2);
for i=1:nCh
    x=[fluorCh{i} '_x'];  y=[fluorCh{i} '_y'];  %metadata x/y field names

    if all( isfield(data.traceMetadata,{x,y}) ),
        coord(i,:) = [data.traceMetadata(m).(x) data.traceMetadata(m).(y)];
    end
end

% Draw markers on selection points.
delete( findall(axFOV,'type','Line') );
viscircles( axFOV, coord, repmat(3,nCh,1), 'EdgeColor','w' );


end



function [p,f,e] = fileparts2(fname)

pidx = find( fname=='\'|fname=='/', 1, 'last' );
eidx = find( fname=='.', 1, 'last' );

p = fname(1:pidx);
f = fname(pidx+1:eidx-1);
e = fname(eidx:end);

end


