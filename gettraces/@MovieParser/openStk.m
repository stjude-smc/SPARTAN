function this = openStk(this, input, params)
%openStk  gettraces internal function to open a movie file.
%
%   STKDATA = openStk(FILENAME, PARMAS)

narginchk(3,3); %FIXME: get params from cascadeConstants if not given.

% Process input arguments, loading a Movie object
if isa(input,'Movie_TIFF') || isa(input,'Movie_STK')
    movie = input;
elseif iscell(input) || ischar(input)
    movie = Movie.load(input);
else
    error('Invalid input: must be filename or Movie object');
end
this.params = params;
this.movie = movie;

% Average the first 10 frames to create an image for finding molecules.
nFrames = movie.nFrames;
averagesize = min([10 nFrames]);

[fields,this.fnames] = subfield( this.movie, params.geometry, 1:averagesize );
fields = cellfun( @(x)mean(x,3), fields, 'Uniform',false );

% Create an estimated background image by sampling the lowest 15-20% of 
% values in each 6x6 area in the image and interpolating values in between.
% The block size (den) should be at least 3x the PSF width.
den = this.params.bgBlurSize;
this.background = cell( size(fields) );
szField = size(fields{1});
temp = zeros( floor(szField/den) );
win = 1:den;
partition = round( 0.167*(den^2) );  %index of sorted pixel to use 

for f=1:numel(fields)
    % Divide image into den-x-den squares and find 16% lowest value in each.
    for i=1:size(temp,1),
        for j=1:size(temp,2),
            sort_temp = fields{f}(den*(i-1)+win, den*(j-1)+win);
            sort_temp = sort( sort_temp(:) );
            temp(i,j) = sort_temp( partition );
        end
    end
    
    % Rescale the image back to the actual size
    this.background{f} = imresize(temp, szField, 'bicubic');
end

% Substract background image
this.stk_top = cellfun( @minus, fields, this.background, 'Uniform',false );

% Use the lowest quartile of intensities from the end of the movie to estimate
% the fundamental level of background noise in each channel.
% This is used in getPeaks for automatically choosing a picking threshold.
endFields = subfield(movie, this.params.geometry, nFrames-11:nFrames-1);
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
[this.total_t, this.peaks, this.total_peaks, this.rejectedTotalPicks,...
 this.fractionOverlapped, this.alignStatus, this.regionIdx, this.psfWidth, ...
 this.integrationEfficiency, this.fractionWinOverlap, this.bgMask] = deal([]);
             
end %FUNCTION OpenStk



