classdef Movie_PMA < Movie
% Movie_PMA  class for loading .pma raw movie files.
%
%   See also: Movie, Movie_STK, Movie_TIFF, openStk, gettraces.

%   Copyright 2018 Cornell University All Rights Reserved.

% The readFrames() way of getting data is a reasonable interface, but it would
% be even cleaner if the data were wrapped almost like mmapPassthrough so that
% it appears to be a data "property" of the movie, but instead just contains
% pointers and methods to read the data. Could call it freadPassthrough. On the
% other hand, am I making things more complicated than needed?

% FIXME: we assume data type is uint16


methods
    
    function obj = Movie_PMA( filename, datatype )
        
        % Get data type and byte size, typically uint16 (2 bytes per pixel)
        if nargin<2 || isempty(datatype)
            options = {'uint8','uint16'};
            datatype = listdlg('PromptString','Select movie data type:', ...
                               'SelectionMode','single', 'ListString',options );
            if isempty(datatype), return; end
            datatype = options{datatype};
        end
        temp = zeros(1,1,datatype); %#ok<PREALL>
        temp = whos('temp');
        datasize = temp.bytes;
        
        % Input may be a cell array, but it must only contain one file.
        % For Movie_TIFF, a file list is possible.
        if iscell(filename)
            assert( numel(filename)==1, 'File lists not allowed!' );
            filename = filename{1};
        end
        
        % Read file header, which only contains frame size in pixels
        obj.filename = filename;
        fid = fopen(filename,'r');
        obj.nX = fread(fid, 1, 'uint16');
        obj.nY = fread(fid, 1, 'uint16');
        fclose(fid);
        
        % Infer number of frames from file size
        frameSize = datasize * obj.nX * obj.nY;  %bytes per frame
        temp = dir(filename);
        approxFrames = (temp.bytes - 2*datasize)/frameSize;
        obj.nFrames = round(approxFrames);
        assert( approxFrames-obj.nFrames<1e-3, 'File corrupted?');
        
        % Use frame numbers since .pma files have no time axis metadata
        obj.timeAxis = 1:obj.nFrames;
        obj.offsets = 2*datasize + frameSize*(0:obj.nFrames-1);
        obj.precision = datatype;
        
    end %constructor
    
    
    % Data access methods. Data are only loaded when needed by these functions.
    function data = readFrames( obj, idx )
        
        % Preallocate space for the output data.
        data = zeros( obj.nY, obj.nX, numel(idx), obj.precision );
        
        fid = fopen( obj.filename, 'r' );
        
        for i=1:numel(idx)
            % Move to start of frame and read data.
            % File is row major. Transpose converts to Matlab-style column
            % major orientation.
            fseek( fid, obj.offsets( idx(i) ), 'bof' );
            data(:,:,i) = fread( fid, [obj.nX obj.nY], ['*' obj.precision] )';
        end
        
        fclose(fid);
    end
    
    function data = readFrame( obj, idx )
        data = obj.readFrames(idx);
    end
    
    
end %public methods



end %class Movie_STK
