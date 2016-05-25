function stkData = OpenStk(filename, params)
%openStk  gettraces internal function to open a movie file.
%
%   STKDATA = openStk(FILENAME, PARMAS)

% Load movie data from file.
movie = Movie.load(filename);
nFrames = movie.nFrames;

% Average the first 10 frames of the movie to create an image (stk_top)
% to use for finding molecules.
averagesize = min([10 nFrames]);
stk_top = movie.readFrames(1:averagesize);
stk_top = mean(stk_top,3);

% Create an estimated background image by:
% 1. Divide the image into den*den squares
% 2. For each square, find the fifth lowest number
% 3. Rescaling these values back to the original image size
% 
% *** den should be much larger than the PSF size. A good rule of thumb is
% at least 3x the size of the width (not std) of the PSF. This is difficult
% to generalize just from the size of the frame.
den = 6;
params.bgBlurSize = den;

background = stk_top;  %**
temp = zeros( floor(movie.nY/den), floor(movie.nX/den) );

for i=1:size(temp,1),
    for j=1:size(temp,2),
        sort_temp = background(den*(i-1)+1:den*i,den*(j-1)+1:den*j);
        sort_temp = sort( sort_temp(:) );

        temp(i,j) = sort_temp( den );  % get the 1/den % smallest value
    end
end

% Rescale the image back to the actual size
background=imresize(temp,[movie.nY movie.nX],'bicubic');

% Background image from the last few frames, used for picking threshold.
endBackground = double( movie.readFrames(nFrames-11:nFrames-1) );
fields = subfield(endBackground,params.geometry);
endBackground = sum( cat(4,fields{:}), 4 );

% Combine the image stack data with the extras just calculated for later
% processing and return them.
% All of this could be combined into the Movie class! TODO
stkData = struct('nChannels',numel(fields), 'movie',movie, 'stk_top',stk_top, ...
                 'background',background,'endBackground',endBackground);

end %FUNCTION OpenStk