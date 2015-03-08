classdef Movie
% This is an abstract class that represents various formats of image stacks
% (movies) recorded from CCD cameras. Format-specific classes inherit from
% this class, creating a standard interface regardless of file format.
% All access is read-only and the file must exist as long as data access is
% needed. Data are not (necessarily) loaded into memory.
%

% Assumptions: the image stack is a series of images with constant exposure
% time and no breaks in imaging?
%
% TO DO: need to define at least some metadata fields that are standardized.
% conversions can be made in the inheriting classes.


properties (SetAccess=protected, GetAccess=public)
    filename; % full path and filename to loaded file - empty if not loaded.
    
    nX;       % size (in pixels) of x dimension (columns).
    nY;       % size (in pixels) of y dimension (rows).
    nFrames;  % number of images in the stack.
    
    timeAxis; % wall time of the start of each frame (starts with zero).
    
end %end public properties


methods (Abstract)
    
    % Data access methods. Data are only loaded when needed by these functions.
    data = readFrames( obj, idxStart, idxEnd );
    data = readFrame(  obj, idx );
    
    
end %public methods



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

        % Load the movie with the appropriate subclass
        [~,~,ext] = fileparts(filename);
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
