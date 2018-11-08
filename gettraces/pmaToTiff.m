function pmaToTiff(filename, precision)
% pmaToTiff  Convert .pma raw frame data to TIFF file image stack
%
%   pmaToTiff(FILENAME)
%   Load a .pma movie file and save a corresponding TIFF-format file
%   compatible with gettraces. If no input is given, the user will be
%   prompted for the location of a file.
%
%   pmaToTiff(FILENAME, PRECISION) specifies the data type of the pixel
%   data, which is typically uint16 (default) or uint8.
%
% See also: Movie, Movie_PMA

%   Copyright 2018 Cornell University All Rights Reserved.


%% Request file path from user if not given as an argument
if nargin<1
    prompt = 'Select a .pma movie file';
    filter = {'*.pma','PMA movie files (*.pma)'; ...
              '*.*','All Files (*.*)'};
    [f,p] = uigetfile(filter, prompt, 'MultiSelect','off');
    
    if f==0, return; end  %user hit cancel
    filename = fullfile(p,f);
end

[p,f] = fileparts(filename);
outFilename = fullfile(p, [f '.tif']);
if exist(outFilename,'file')
    delete(outFilename);
end

if nargin<2, precision = 'uint8'; end


%% Load PMA input file
% FIXME: use Tiff class if available for faster writing?

fid = fopen(filename,'r');
nX = fread(fid, 1, 'uint16');  %frame size in pixels
nY = fread(fid, 1, 'uint16'); 

% Load frames
while ~feof(fid)
    frame = fread(fid, [nX nY], ['*' precision])';
    if isempty(frame), break; end
    imwrite(frame, outFilename, 'WriteMode', 'append',  'Compression','none' );
end


end %function pmaToTiff

