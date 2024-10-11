function this = openStk(this, input, params)
% MovieParser function to open a new movie
% input can be a Movie object or the path to a movie file.
% params comes from cascadeConstants.m.

narginchk(2,3);

if nargin>=3
    this.params = params;
end

% Load movie from file
if isa(input,'Movie')
    movie = input;
elseif iscell(input) || ischar(input)
    movie = Movie.load(input);
else
    error('Invalid input: must be filename or Movie object');
end

% NOTE: the movie metadata may include significantly more information; we
% only use the channel names to keep it simple for now.
metadataFound = false;

% If available and valid, use movie's metadata to assign channels.
if this.autoFromMetadata && ...
   all( isfield(movie.metadata,{'fieldArrangement','channels'}) ) && ...
   all( isfield(movie.metadata.channels,{'name','wavelength'})  )
    
    profileNames = {this.params.channels.name}; %all possible channel names for loaded profile.
    metaNames = {movie.metadata.channels.name}; %names listed in movie metadata
    geo = movie.metadata.fieldArrangement;
    
    if ~all( ismember(metaNames,profileNames) )
        warning('Movie metadata lists channel names that do not match loaded profile! Ignoring');
    elseif ~issorted([movie.metadata.channels.wavelength])
        warning('Movie metadata channels not in wavelength order! Ignoring');
    elseif ~all( geo(:)>=0 & geo(:)<=numel(metaNames) & geo(:)==floor(geo(:)) & sum(geo(:)>0)==numel(metaNames) )
        warning('Movie metadata field arrangement invalid! Ignoring');
    else
        metadataFound = true;
        
        % Create object for extracting spectral channels from frame data.
        % photonsPerCount from metadata overrides profile.
        [~,this.idxActiveChannels] = ismember(metaNames, profileNames);
        %ch = this.params.channels(this.idxActiveChannels);
        ch = movie.metadata.channels;
        this.chExtractor = ChannelExtractor(movie, geo, ch);

        % If no roles previously assigned or the number of channels has
        % changed, make a reasonable guess. This will also be done in
        % whatever function allows the user to change the geometry.
        nCh = numel( this.idxActiveChannels );
        if numel(this.roles)~=nCh
            this.roles = guessRoles(nCh);
        end
    end
end %if metadata found
    
% If metadata is not available, use existing settings to load new movie.
% FIXME: this may give an error if the movie isn't evenly divisible.
if ~metadataFound && ~isempty(this.chExtractor)
    this.chExtractor.movie = movie;

% If this is the first run and no metadata available, guess.
elseif isempty(this.chExtractor)
    geo = guessFieldArrangement(movie, this.params.geometry, this.params.channels);
    this.idxActiveChannels = sort( geo(geo>0) );
    ch = this.params.channels(this.idxActiveChannels);
    
    for i=1:numel(ch)
        geo( geo==this.idxActiveChannels(i) ) = i;
    end
    this.chExtractor = ChannelExtractor(movie, geo, ch);
    this.roles = guessRoles( numel(ch) );
end

% Average the first few frames to create an image for finding molecules.
this.chExtractor.avgTop( this.params.nAvgFrames, this.params.subtractBGImage );

% Use the lowest quartile of intensities from the end of the movie to estimate
% the fundamental level of background noise in each channel.
% This is used in getPeaks for automatically choosing a picking threshold.
% FIXME: in future, get std from each field separately, otherwise difference in
% the background level (bias) dominate. SEE NEXT SECTION
s = max(1, this.nFrames-11);
endFields = this.chExtractor.read( s:this.nFrames-1 );
endBG = sum( cat(4,endFields{:}), 4 );
endBG = sort(endBG(:));
endBG = endBG( 1:floor(numel(endBG)*0.75) );
this.stdbg = std(endBG);

% Improved version that requires a different threshold
% this.stdbg = zeros( numel(fields),1 );
% for f=1:numel(endFields)  %better version
%     sort_temp = sort( endFields{f}(:) );
%     sort_temp = sort_temp( 1:floor(numel(sort_temp)*0.75) );
%     this.stdbg(f) = std( double(sort_temp) );
% end
        
% Reset any stale data from later steps
[this.total_t, this.peaks, this.total_peaks, this.rejectedPicks, this.rejectedTotalPicks,...
 this.fractionOverlapped, this.alignStatus, this.regionIdx, this.bgMask, ...
 this.integrationEfficiency, this.psfWidth, this.fractionWinOverlap ] = deal([]);


end %FUNCTION OpenStk


function roles = guessRoles(nCh)
% Make an educated guess of channel roles just based on the number.
    if nCh<4
        temp = {'donor','acceptor','acceptor2'};
        roles = temp(1:nCh);
    else
        roles = {'donor','acceptor','donor2','acceptor2'};
    end
end
