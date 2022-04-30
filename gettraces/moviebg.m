function background = moviebg(fields, bgBlurSize)

narginchk(1,2);
nargoutchk(1,1);

if isnumeric(fields),
    fields = {fields};
end
if ~iscell(fields)
    error('Invalid input type');
end

if nargin<2
    %params = cascadeConstants();
    bgBlurSize = 6;  %params.bgBlurSize;
end

% Average all frames input
fields = cellfun( @(x)mean(x,3), fields, 'Uniform',false );

% Create an estimated background image by sampling the lowest 15-20% of 
% values in each 6x6 area in the image and interpolating values in between.
% The block size (den) should be at least 3x the PSF width.
den = bgBlurSize;
background = cell( size(fields) );
szField = size(fields{1});
temp = zeros( floor(szField/den) );
win = 1:den;
partition = floor( 0.167*(den^2) );  %index of sorted pixel to use

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
    background{f} = imresize(temp, szField, 'bicubic');
end

end %function
