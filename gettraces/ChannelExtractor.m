classdef ChannelExtractor < handle
% Demultiplex fluorescence spectral channels from raw camera frame data.
%
%   Raw camera-based smFRET data is recorded as a set of images (frames),
%   each of which is a 2-dimension matrix of pixel intensities.
%   These images often contain multiple projections of the field of view,
%   each corresponding to a distinct fluorophore's emission (spectral band).
%   This class provides an interface to convert raw frame data into these
%   spectral "channels".
%
%   chExtr = ChannelExtractor(movie) creates a new ChannelExtractor chExtr,
%   where movie is a path to a movie file or a Movie object. If the movie's
%   metadata encodes channel information, this will be used to split the
%   pixel data into channels. If not, assume single channel (Cy3).
%
%   chExtr = ChannelExtractor( movie, fieldArranagement, channels ) 
%   specify all information needed to split movie data into channels.
%    * fieldArranagement (int): matrix whose shape describes how to divide frame 
%       data into spectral channels; elements specify the channel order.
%    * channels: struct array (one per channel) with fields:
%       - name (cell array of strings): simple name for each channel.
%       - wavelength (int): laser line for each channel (for display only).
%       - photonsPerCount (double): camera units to photons conversion (optional).
%       - role (string): donor, acceptor, ignore, etc. (optional).
%
%   Examples:
%      ch = struct('name',{'Cy3','Cy5'}, 'wavelength',{532,640}, 'photonsPerCount',0.49);
%      chExtr = ChannelExtractor( 'E:\movieLR.tif',  [1 2],      ch );   % 2-color (L,R)
%      chExtr = ChannelExtractor( 'E:\movieCAT.tif', cat(3,1,2), ch );   % 2-color (concatinated frames)
%
%   channelCellArray = chExtr.read(1:10) reads the first 10 frames of the
%   movie and returns a cell array, one element per channel, each
%   is a matrix with image data (dimensions: columns, rows, frames).
%
%   NOTE: all properties are public, but if changing the number of
%   channels, it is best to call the "verify" method.
%
% See also: Movie, TraceExtractor, gettraces_gui

% Copyright 2022 All Rights Reserved.



%%
% Any correction parameters (scaling, crosstalk, ade) provided here should
% be seen as defaults for initial analysis to be refined in later steps.

properties (SetAccess=public, GetAccess=public)

    movie;          % Movie object for reading raw frame data from file.
    fieldArranagement = 1;   % matrix describing how to split frame data into channels.
    
    % struct array describing spectral channels (N):
    % - name (string): generic identifier of the spectral band (e.g., Cy3).
    % - description (string): detail about the sample being imaged.
    % - wavelength (int): emission max or relevant laser line in nm.
    % - photonsPerCount (dbl): conversion from camera units to detected photons.
    % - role (string): double, acceptor, ignore, etc.
    channels = struct( 'name','Cy3', 'description','', 'wavelength',532, ...
                       'photonsPerCount',1.0, 'role','' );
    
    % Non-standard settings
    skipFrames = 0;   %number of frames to ignore at the beginning.
    %cropX;   %remove left and right pixels so image is cropX cols wide.
    %cropY;   %remove top and bottom pixels so image is cropX rows tall.
    %interleaved;  %if true, channels are interleaved over frames.
end

% These are calculated from the input parameters
properties (Dependent)
    nFrames;
    nChannels;
    nX;
    nY;
    timeAxis;
end


%%
methods
    
    % Constructor
    function this = ChannelExtractor(varargin)
        this.load(varargin{:});
    end
    
    
    function load(this, varargin)
        narginchk(2,4);
        required = {'name','wavelength'};
        
        % First argument: movie file path or Movie object.
        input = varargin{1};
        if ischar(input) || iscell(input)
            this.movie = Movie.load(input);
        elseif isa(input,'Movie') && numel(input)==1
            this.movie = input;
        else
            error('Invalid first argument. Must be file path or Movie object');
        end
        
        % Load metadata from movie file if available.
        % If no metadata, assume single-channel (Cy3 donor etc).
        if nargin<1 && isfield(this.movie.metadata,'channels')
            initCh = this.channels;
            try
                this.fieldArranagement = this.movie.metadata.fieldArrangement;
                this.channels = mergeStruct( this.movie.metadata.channels, this.channels, required );
            catch
                warning('Movie metadata found but invalid. Defaulting to single channel.');
                this.channels  = initCh;
                this.fieldArranagement  = 1;
            end
        end
        
        % Parse field arrangement and channel list parameters
        if nargin>2
            narginchk(4,4);
            assert( isnumeric(varargin{2}), 'Invalid second input. Should be a numeric matrix' );
            this.fieldArranagement = varargin{2};

            assert( isstruct(varargin{3}), 'Invalid third input. Should be a struct array' );
            this.channels = mergeStruct( varargin{3}, this.channels, required );
        end
        
        % Verify all properties are self-consistent.
        this.verify();
        
    end %function load
    
    
    
    function output = read( this, idx )
        % This function returns a cell array. Each element corresponds to
        % one camera (spectral channel) and contains the frame data as a 3D
        % array (rows, columns, frames).
        
        idx = idx + this.skipFrames;
        assert( all(idx>=1 & idx<=this.nFrames & idx==floor(idx)), 'Invalid index' );

        % Parse image data into channels
        output = splitFrame( this.movie, this.fieldArranagement, idx );
        
    end %function readFrames
    
    
    
    % Check internal state and give an error if invalid
    function verify(this)
        assert( isa(this.movie,'Movie') && numel(this.movie)==1 );
        
        % Verify field arrangement make sense with channels struct.
        szCh = size(this.fieldArranagement);
        if numel(szCh)>2 && szCh(1:2)>1
            error('Channels may be tiled in space or as concatinated frames, but not both.');
        end
        assert( all(ismember(this.fieldArranagement(:),0:this.nChannels)), 'Invalid fieldArranagement' );
        ch = this.fieldArranagement(this.fieldArranagement>0);
        assert( numel(unique(ch))==this.nChannels, 'Duplicate channels not allowed' )
        
        % Verify channel fields are single values with correct type.
        for i=1:numel(this.channels)
            ch = this.channels(i);
            assert( ischar(ch.name) );
            assert( numel(ch.wavelength)==1 && isnumeric(ch.wavelength) );
            assert( numel(ch.photonsPerCount)==1 && isnumeric(ch.photonsPerCount) );
        end
    end
    
    
    
    % Number of actual movie frames, after splitting files that have
    % movies from multiple cameras concatinated together.
    function output = get.nFrames(this)
        output = (this.movie.nFrames-this.skipFrames) / size(this.fieldArranagement,3);
    end
    
    function output = get.nChannels(this)
        output = numel(this.channels);
    end
    
    function output = get.nX(this)
        output = this.movie.nX;
    end
    
    function output = get.nY(this)
        output = this.movie.nY;
    end
    
    function output = get.timeAxis(this)
        %NOTE: assume concatinated (not interleaved) frames
        output = this.movie.timeAxis( (1+this.skipFrames):this.nFrames );
    end
    
end  %methods




end %class




function template = mergeStruct( input, template, required )
% Copy data from INPUT into TEMPLATE, ignoring any fields not in TEMPLATE.

    % Expand template if input is a struct array.
    if numel(input)>1
        assert( numel(template)==1 );
        template = repmat( template, size(input) );
    end

    % Copy data from input struct, ignoring extraneous fields.
    fn = fieldnames(input);
    fn = fn( ismember(fn,fieldnames(template)) );
    for i=1:numel(fn)
        [template.(fn{i})] = input.(fn{i});
    end
    
    % Verify all required fields have been provided
    if nargin>=3 && ~all(ismember(required,fn))
        error('Not all required fields given for struct input');
    end

end

