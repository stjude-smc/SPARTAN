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


% NOTES
% - Rename class to ChannelManager and read() to readFrames.
% - This class will store and manage the current profile (as it relates to
%   channels) and provide methods for updating values, including correction parameters.
% - set movie method will merge current state with the new metadata and
%   give an error if it doesn't match.
% - 


%%
% Any correction parameters (scaling, crosstalk, ade) provided here should
% be seen as defaults for initial analysis to be refined in later steps.

properties (SetAccess=public, GetAccess=public)

    movie;                  % Movie object for reading raw frame data from file.
    fieldArrangement = 1;   % matrix describing how to split frame data into channels.
    
    % struct array describing spectral channels (N):
    % - name (string): generic identifier of the spectral band (e.g., Cy3).
    % - description (string): detail about the sample being imaged.
    % - wavelength (int): emission max or relevant laser line in nm.
    % - photonsPerCount (dbl): conversion from camera units to detected photons.
    channels = struct( 'name','Cy3', 'description','', 'wavelength',532, ...
                       'photonsPerCount',1.0 );
                   
    % struct array describing illumination periods (ALEX)
    lasers = [];
    
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
        
        % Default value for illumination metadata.
        %this.lasers = struct('wavelength',532, 'framesActive',1:this.movie.nFrames);
        
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
        if isfield(this.movie.metadata,'lasers') && ...
           all( isfield(this.movie.metadata.lasers,{'wavelength','framesActive'}) )
            assert( issorted([this.movie.metadata.lasers.wavelength]) );
            this.lasers = this.movie.metadata.lasers;
        end
        
        if nargin<2 && all( isfield(this.movie.metadata,{'channels','fieldArrangement'}) )
            initCh = this.channels;
            try
                this.fieldArrangement = this.movie.metadata.fieldArrangement;
                this.channels = mergeStruct( this.movie.metadata.channels, this.channels, required );
            catch
                warning('Movie metadata found but invalid. Defaulting to single channel.');
                this.channels = initCh;
                this.fieldArranagement = 1;
            end
        
        % Parse field arrangement and channel list parameters
        elseif nargin>2
            narginchk(4,5);
            assert( isnumeric(varargin{2}), 'Invalid second input. Should be a numeric matrix' );
            this.fieldArrangement = varargin{2};

            assert( isstruct(varargin{3}), 'Invalid third input. Should be a struct array' );
            this.channels = mergeStruct( varargin{3}, this.channels, required );
            
            if nargin>=5
                this.lasers = varargin{4};
            end
        end
        
        % Verify all properties are self-consistent.
        this.verify();
        
    end %function load
    
    
    function wavelengths = lasersActive(this, idx)
        % Return which laser wavelengths are active in frame IDX
        wavelengths = [];
        if isempty(this.lasers), return; end
        
        for i=1:numel(this.lasers)
            if any(this.lasers(i).framesActive==idx)
                wavelengths(end+1) = this.lasers(i).wavelength; %#ok<AGROW>
            end
        end
    end
    
    
    function output = read( this, idx )
        % This function returns a cell array. Each element corresponds to
        % one camera (spectral channel) and contains the frame data as a 3D
        % array (rows, columns, frames).
        
        idx = idx + this.skipFrames;
        assert( all(idx>=1 & idx<=this.nFrames & idx==floor(idx)), 'Invalid index' );

        % Parse image data into channels
        output = splitFrame( this.movie, this.fieldArrangement, idx );
        
    end %function readFrames
    
    
    
    % Check internal state and give an error if invalid
    function verify(this)
        assert( isa(this.movie,'Movie') && numel(this.movie)==1 );
        
        % Verify field arrangement make sense with channels struct.
        szCh = size(this.fieldArrangement);
        if numel(szCh)>2 && szCh(1:2)>1
            error('Channels may be tiled in space or as concatinated frames, but not both.');
        end
        assert( all(ismember(this.fieldArrangement(:),0:this.nChannels)), 'Invalid fieldArrangement' );
        ch = this.fieldArrangement(this.fieldArrangement>0);
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
        output = (this.movie.nFrames-this.skipFrames) / size(this.fieldArrangement,3);
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

