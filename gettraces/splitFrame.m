function output = splitFrame(input, geo, frameIdx)
%splitFrame   Split image data into uniform-size subfields
%
%   OUT = splitFrame(MOVIE,GEO,FRAMES) loads the frame numbers specified in
%   FRAMES and divides these images into equal-sized subregions where the 
%   shape of the array GEO identifies along which dimensions to divide the
%   images. The binary values in GEO define which subregions to return in
%   the cell array OUT.
%   
%   GEO examples:
%     [1 0; 0 2] split image 2x2 and select upper-left and lower-right.
%     cat(3,1,2) images stacked sequentially (2 channels).

%   Copyright 2016-2022 All Rights Reserved.


% Process input arguments
narginchk(3,3);
assert( isa(input,'Movie'), 'Invalid input' );
assert( numel(size(geo))<=3, 'Invalid dimensions of subfield indexing matrix' );

% Fluorescence channels are stacked in sequential frames (not interleaved)
if size(geo,3)>1
    if size(geo,1)>1 || size(geo,2)>1
        error('Fields must be arranged either side-by-side OR stacked, not both');
    end
    
    nFrames = input.nFrames/size(geo,3);  %number of actual time units in movie
    assert( nFrames==floor(nFrames), 'Unexpected number of frames for chosen stacked geometry' )
    
    output = cell( sum(geo>0), 1 );
    
    for i=1:size(geo,3)
        if geo(i)==0, continue; end
        output{ geo(i) } = input.readFrames( frameIdx + (i-1)*nFrames );
    end
    
% Fluorescence channels are stiched side-by-side within each frame
else
    input = input.readFrames(frameIdx);

    % Determine number and size of subfields
    [nr,nc] = size(geo);
    [imr,imc,imf] = size(input);
    idx = [imr,imc] ./ [nr,nc];  %subdivision size in each dimension

    if ~all(idx==floor(idx))
        error('Movie cannot be divided into equal-sized fields. Incorrect geometry?');
    end

    % Divide image into equal-sized subregions
    if imf==1
        C = mat2cell( input, repmat(idx(1),nr,1), repmat(idx(2),nc,1) );
    else
        C = mat2cell( input, repmat(idx(1),nr,1), repmat(idx(2),nc,1), imf );
    end
    
    output = cell( sum(geo(:)>0), 1 );
    for i=1:numel(geo)
        if geo(i)>0
            output( geo(i) ) = C(i);
        end
    end
end

end %FUNCTION subfield

