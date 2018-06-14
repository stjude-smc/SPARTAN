classdef Movie_STK < Movie
% Movie_STK is a wrapper for a stack of TIFF images, as generated by imaging
% software like MetaMorph. See Movie for interface details.

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% TODO: first get this class to work with loading headers with every frame
% having a seperate entry. Then get everything working with gettraces. Then
% try to optimize the code by, for example, only reading the first frame's
% annotation.
%
% The readFrames() way of getting data is a reasonable interface, but it would
% be even cleaner if the data were wrapped almost like mmapPassthrough so that
% it appears to be a data "property" of the movie, but instead just contains
% pointers and methods to read the data. Could call it freadPassthrough. On the
% other hand, am I making things more complicated than needed?


properties (SetAccess=protected, GetAccess=public)
    % See Movie class for additional standard properties.
end %end public properties


properties (SetAccess=protected, GetAccess=protected),
    dataOffsets = []; %offsets to data segments in the movie file, with
                      %one offset per plane (frame).
end


methods
    
    function obj = Movie_STK( filename )
        
        % Input may be a cell array, but it must only contain one file.
        % For Movie_TIFF, a file list is possible.
        if iscell(filename)
            assert( numel(filename)==1, 'Movie_STK: file lists not allowed!' );
            filename = filename{1};
        end
        
        obj.filename = filename;
        
        % Read the TIFF file and get all useful header information. The data
        % sections are ignored. Data-access offsets are stored in tiffData.
        [obj.header,obj.dataOffsets] = readTiffHeader( filename );
        
        % Extract basic image metadata
        obj.nX = obj.header(1).width;
        obj.nY = obj.header(1).height;
        obj.nFrames = numel( obj.dataOffsets );
        
        % Generate an approximate time axis. The actual timestamps are in the
        % MM_private1 (UIC1, 33628) field untag tag #16 (CreateTime). The LONG
        % data element is a pointer to a LONG [date,time]. FIXME
        x = repmat( obj.header(1).MM.Exposure, [1 obj.nFrames] );
        obj.timeAxis = [0 cumsum(x(1:end-1))];
        
    end %constructor
    
    
    % Data access methods. Data are only loaded when needed by these functions.
    function data = readFrames( obj, idx )
        
        % Parse input arguments; insure correct orientation of vector.
        idx = reshape(idx, [1 numel(idx)]);
        assert( min(idx)>=1 && max(idx)<=obj.nFrames, 'Invalid indexes' );
        
        % Preallocate space for the output data.
        data = zeros( obj.nY,obj.nX,numel(idx), 'uint16' );
                
        fid = fopen( obj.filename );
        
        framesRead=0;
        for i=idx,
            % Move to the start of the data section
            fseek( fid, obj.dataOffsets(i), 'bof' );
            
            % FIXME: supports only uint16!!
            frame = fread( fid, [obj.nX obj.nY], '*uint16' );
            data(:,:,framesRead+1) = frame';
            
            framesRead=framesRead+1;
        end
        
        fclose(fid);
    end
    
    function data = readFrame( obj, idx )
        fid = fopen( obj.filename );
        fseek( fid, obj.dataOffsets(idx), 'bof' );
        data = fread( fid, [obj.nX obj.nY], '*uint16' )';
        fclose(fid);
    end
    
    
end %public methods



end %class Movie_STK
