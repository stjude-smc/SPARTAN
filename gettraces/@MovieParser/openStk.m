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

% FIXME: need to re-evaluate the paradigm here. This class ends up doing a
% ton of manipulation of ChannelExtractor's properties. Does that suggest
% that this code belongs in that class somehow? Or does it suggest this
% property should be in this class instead?

% If available and valid, use movie's metadata to assign channels.
% FIXME: assign any fields not in metadata from params.
if all( isfield(movie.metadata,{'fieldArrangement','channels'}) )
    geo = movie.metadata.fieldArrangement;
    ch = movie.metadata.channels;
    
    if all(ismember({ch.name},{this.params.channels.name}))
        % Replace current channel settings with those from file.
        try
            this.chExtractor = ChannelExtractor(movie, geo, ch);
        catch e
            warning(e.identifier,'Failed to load channel metadata from movie: %s',e.message);
            this.chExtractor = [];
        end
    else
       warning('Movie metadata and current profile not compatible');
       this.chExtractor = [];
    end

% If no metadata is available in the file, but settings were defined for
% the previous movie, try to retain them if possible.
elseif ~isempty(this.chExtractor)
    geo = this.chExtractor.fieldArranagement;
    ch  = this.chExtractor.channels;
    this.chExtractor = ChannelExtractor(movie,geo,ch);
end

% If this is the first run or metadata was invalid, default to single channel.
if isempty(this.chExtractor)
    geo = 1;
    [~,idx] = min(  abs([this.params.channels.wavelength]-532)  );  %arbitrary
    ch = this.params.channels(idx);
    this.chExtractor = ChannelExtractor(movie, geo, ch);
end

% If no roles given/retained, provide a reasonable guess.
if ~isfield(this.chExtractor.channels,'role') || isempty(this.chExtractor.channels(1).role)
    if this.chExtractor.nChannels<4
        roles = {'donor','acceptor','acceptor2'};
    else
        roles = {'donor','acceptor','donor2','acceptor2'};
    end
    [this.chExtractor.channels.role] = roles{ 1:this.chExtractor.nChannels };
end

% Average the first few frames to create an image for finding molecules.
averagesize = min( [this.params.nAvgFrames this.nFrames] );
fields = this.chExtractor.read( 1:averagesize );
fields = cellfun( @(x)mean(x,3), fields, 'Uniform',false );

% Substract background image
this.background = moviebg(fields);
this.stk_top = cellfun( @minus, fields, this.background, 'Uniform',false );

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



