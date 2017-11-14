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

% Average the first 10 frames of the movie to create an image (stk_top)
% to use for finding molecules.
nFrames = movie.nFrames;
averagesize = min([10 nFrames]);
% this.stk_top = mean( movie.readFrames(1:averagesize), 3);

% Extract individual fluorescence fields
this.fields = subfield( this.movie, params.geometry, 1:averagesize );
this.fields = cellfun( @(x)mean(x,3), this.fields, 'Uniform',false );
this.nChannels = numel(this.fields);

% Background subtracted version 
this.stk_top = cellfun( @minus, stkData.fields, stkData.background, 'Uniform',false );

% Create an estimated background image by:
% 1. Divide the image into den*den squares
% 2. For each square, find the fifth lowest number
% 3. Rescaling these values back to the original image size
% 
% *** den should be much larger than the PSF size. A good rule of thumb is
% at least 3x the size of the width (not std) of the PSF. This is difficult
% to generalize just from the size of the frame.
den = 6;
this.params.bgBlurSize = den;  %FIXME
this.background = cell( size(this.fields) );
szField = size(this.fields{1});
temp = zeros( floor(szField/den) );  %movie.nY,nX

for f=1:numel(this.fields)
    
    for i=1:size(temp,1),
        for j=1:size(temp,2),
            sort_temp = this.fields{f}(den*(i-1)+1:den*i,den*(j-1)+1:den*j);
            sort_temp = sort( sort_temp(:) );

            temp(i,j) = sort_temp( den );  %get the 1/den % smallest value
        end
    end
    
    % Rescale the image back to the actual size
    this.background{f} = imresize(temp, szField, 'bicubic');
end

% Background image from the last few frames, used for picking threshold.
endFields = subfield(movie, this.params.geometry, nFrames-11:nFrames-1);
this.endBackground = sum( cat(4,endFields{:}), 4 );  %fixme: cell array instead??
        
% Reset any stale data from later steps
this.stage = 1;
[this.total_t, this.peaks, this.total_peaks, this.rejectedTotalPicks,...
 this.fractionOverlapped, this.alignStatus, this.regionIdx,...
 this.integrationEfficiency, this.fractionWinOverlap, this.bgMask] = deal([]);
             
end %FUNCTION OpenStk



