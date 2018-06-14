classdef Movie < handle
% This is an abstract class that represents various formats of image stacks
% (movies) recorded from CCD cameras. Format-specific classes inherit from
% this class, creating a standard interface regardless of file format.
% All access is read-only and the file must exist as long as data access is
% needed. Data are not (necessarily) loaded into memory.
%
% See also: Movie_STK, Movie_TIFF, gettraces.

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% Assumptions: the image stack is a series of images with constant exposure
% time and no breaks in imaging?
%
% TO DO: need to define at least some metadata fields that are standardized.
% conversions can be made in the inheriting classes.


properties (SetAccess=protected, GetAccess=public)
    filename; % full path and filename to loaded file - empty if not loaded.
    
    nX=0;       % size (in pixels) of x dimension (columns).
    nY=0;       % size (in pixels) of y dimension (rows).
    nFrames=0;  % number of images in the stack.
    
    timeAxis=[];   % wall time of the start of each frame (starts with zero).
    
    header = struct([]); %
    
end %end public properties


methods (Abstract)
    % Data access methods
    data = readFrames( obj, idxStart, idxEnd );
end


methods
    function viewer = show(this)
    % Open a MovieViewer window for user to display the movie
    viewer = MovieViewer(this);
    viewer.show();
    end
end


methods (Static)
    
    % Movie-making factory method.
    function obj = load( filename )
        % Ask the user for a file if none given.
        if nargin<1 || isempty(filename),
            [f,p] = uigetfile( '*.stk;*.tif;*.tiff', 'Choose a movie file' );
            if f==0,
                obj = []; %user hit cancel
                return;
            end
            filename = fullfile(p,f);
        end

        if ~iscell(filename),
            filename = {filename};
        end

        % Load the movie with the appropriate subclass
        [~,~,ext] = fileparts(filename{1});
        if ~isempty( strfind(ext,'tif') ),
            obj = Movie_TIFF(filename);
        elseif strcmp(ext,'.stk'),
            obj = Movie_STK(filename);
        else
            error('Unrecognized movie type');
        end
    end
    
end



end %class Movie.
