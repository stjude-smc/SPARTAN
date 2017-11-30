function output = subfield(input, quad, frameIdx)
%SUBFIELD  Extract image subfields designated by a string
%
%   OUT = subfield(MOVIE,QUAD,FRAMES) loads the frame numbers specified in
%   FRAMES and divides these images into equal-sized subregions where the shape
%   of the logical array QUAD identifies along which dimensions to divide the
%   images. The binary values in QUAD define which subregions to return in the
%   cell array OUT.
%   
%   EXAMPLES:
%     If QUAD is logical([1 0; 1 0]), images are split in four regions  
%     and the upper-left and lower-left quadrants are returned in OUT.
%   
%     If QUAD is true(1,1,3), there are three fluorescence channels that 
%     are stacked sequentially in the movie, rather than side-by-side.

%   Copyright 2016-2017 Cornell University All Rights Reserved.


% Process input arguments
narginchk(3,3);
assert( isa(input,'Movie_TIFF') || isa(input,'Movie_STK'), 'Invalid input' );

szQuad = size(quad);
assert( numel(szQuad)<=3, 'Invalid dimensions of subfield indexing matrix' );

% Fluorescence channels are stacked in sequential frames (not interleaved)
if numel(szQuad)>2
    if any(szQuad(1:2)>1)
        error('Fields must be arranged either side-by-side OR stacked, not both');
    end
    
    nFrames = input.nFrames/szQuad(3);  %number of actual time units in movie
    assert( nFrames==floor(nFrames), 'Unexpected number of frames for chosen stacked geometry' )
    
    output = cell( sum(quad>0), 1 );
    for i=to_row( find(quad>0) )
        output{i} = input.readFrames( frameIdx + (i-1)*nFrames );
    end
    
% Fluorescence channels are stiched side-by-side within each frame
else
    input = input.readFrames(frameIdx);

    % Determine number and size of subfields
    [nr,nc] = size(quad);
    [imr,imc,imf] = size(input);
    idx = [imr,imc] ./ [nr,nc];  %subdivision size in each dimension

    if ~all(idx==floor(idx))
        error('Movie cannot be divided into equal-sized fields. Incorrect geometry?');
    end

    % Divide image into equal-sized subregions
    C = mat2cell( input, repmat(idx(1),nr,1), repmat(idx(2),nc,1), imf );
    output = C(quad>0);  %Select requested subfields
end

% Sort fields in wavelength order
if ~islogical(quad)
    idxFields = to_row( quad(quad>0) );
    output = output(idxFields);  %sort in wavelength order
end

end %FUNCTION subfield


