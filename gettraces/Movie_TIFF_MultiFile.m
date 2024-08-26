classdef Movie_TIFF_MultiFile < Movie
% Like Movie_TIFF but will take a cell array of filenames and treat them as
% one larger virtual stack.
%
%         movie = Movie_TIFF_MultiFile( filenames )
% 
% The readFrame(idx) and readFrames(indexes) methods can be used to read
% data from the movie file.

%   Copyright 2007-2024 All Rights Reserved.



%% ============ PROPERTIES ============= %
properties (SetAccess=protected, GetAccess=protected)
    submovies;  %array of Movie_TIFF objects representing input files.
    nFramesPerMovie;  
end



%% ============ PUBLIC METHODS ============= %
methods
    % CONSTRUCTOR: input is a cell array paths to TIFF-format files
    % thta are assumed to already be in the correct order.
    function obj = Movie_TIFF_MultiFile( filenames )
        assert(iscell(filenames) && ~isempty(filenames));
        
        % Load TIFF files and sort by movie start time.
        obj.submovies = cellfun( @(x)Movie_TIFF(x), filenames );

        obj.filename = obj.submovies(1).filename;
        obj.header = obj.submovies(1).header;
        obj.nX = obj.submovies(1).nX;
        obj.nY = obj.submovies(1).nY;
        obj.precision = obj.submovies(1).precision;
        obj.metadata  = obj.submovies(1).metadata;
        obj.nFrames = sum([obj.submovies.nFrames]);

        % Construct a continuous time axis from time interval, ignoring any
        % disrepencies in the metadata.
        if obj.submovies(1).timeAxis(1)==1
            obj.timeAxis = 1:obj.nFrames;
        else
            frame_interval_ms = diff(obj.submovies(1).timeAxis(1:2));
            obj.timeAxis = frame_interval_ms*(0:obj.nFrames-1);
        end

    end %constructor
    
    
    function data = readFrames( obj, input )
        %assert( numel(idx)==1 && idx>=1 && idx<=obj.nFrames, 'Invalid index' ); 
        movieFirstFrame = 1+cumsum([0 obj.submovies.nFrames]);
        data = zeros( obj.nY,obj.nX,numel(input), obj.precision );

        for i=1:numel(input)
            idxFile = find( input(i)>=movieFirstFrame, 1, 'last' );
            idx = input(i) - movieFirstFrame(idxFile) +1;
            data(:,:,i) = obj.submovies(idxFile).readFrames(idx);
        end
    end
    
end %public methods



end %class

