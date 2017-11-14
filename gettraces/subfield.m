function output = subfield(input, quad, frameIdx)
%SUBFIELD  Extract image subfields designated by a string
%
%   OUT = subfield(IN,STR) extracts a sub-fields from the image matrix IN
%   as described in STR, which can be empty (whole field), L (left), R (right),
%   T (top), B (bottom), TL, TR, BL, or BR. If STR is a cell array, OUT will
%   be a cell array with a subfield associated with each cell in STR.
%   
%   OUT = subfield(IN,GEO) extracts standard set of sub-fields. GEO can be:
%   1 (whole field), 2 (L/R), 3 (T/B), or 4 (quadrants).

%   Copyright 2016-2017 Cornell University All Rights Reserved.


% Process input arguments
narginchk(3,3);
nargoutchk(1,1);
assert( isa(input,'Movie_TIFF') || isa(input,'Movie_STK'), 'Invalid movie input' );  %why do they not both inheret Movie?
if ~isnumeric(frameIdx) || ~isvector(frameIdx) || any(frameIdx<1|frameIdx>input.nFrames),
    error('Invalid frame index');
end

% Single string targeting a field to extract
if ischar(quad),
    output = subfield2(input,quad,frameIdx);
    
% Cell array of strings targeting each field to extract
elseif iscell(quad)
    output = cell(size(quad));
    for i=1:numel(quad)
        output{i} = subfield2(input,quad{i},frameIdx);
    end
    
    % Truncate fields to the same size (should never happen)
    sz = cellfun(@size, output, 'Uniform',false);
    
    if ~all(  cellfun(@(x)isequal(sz{1},x), sz)  )
        warning('gettraces:fieldSizeMismatch','Movie cannot be divided into equal-sized subfields. Incorrect channel geometry? Check settings.');
        
        szMin = min( cat(1,sz{:}) );
        for i=1:numel(output)
            output{i} = output{i}( 1:szMin(1), 1:szMin(2) );
        end
    end
    
% Integer specifing geometry. Return pre-defined set of sub-images.
elseif isnumeric(quad)
    switch quad
        case 1,  output = { double(input.readFrames(idxFrames)) };
        case 2,  output = subfield(input, {'L','R'}, frameIdx);
        case 3,  output = subfield(input, {'T','B'}, frameIdx);
        case 4,  output = subfield(input, {'TL','TR','BL','BR'}, frameIdx);
    end
else
    error('Invalid field identifier');
end

end %FUNCTION subfield



%%
function output = subfield2(input, quad, frameIdx)
% Extract subfield for scalar input

if isempty(quad)
    output = input.readFrames(frameIdx);

elseif quad(1)=='S'||quad(1)=='I'
    idx = sscanf(quad(2:end),'%d/%d');
    assert(idx(2)>idx(1) & idx(2)<=4 & idx(1)>0);
    nFrames = input.nFrames/idx(2);
    
    % Fields are pasted sequentially in time
    if quad(1)=='S'
        frameIdx = frameIdx + nFrames*(idx(1)-1);

    % Fields are interleaved
    elseif quad(1)=='I'
        chFrames = idx(1):idx(2):input.nFrames;
        frameIdx = chFrames(frameIdx);
    end
    output = input.readFrames(frameIdx); 

% Fields are stitched side-by-side. Each frame has all channels.
else
    input = input.readFrames(frameIdx);
    [nrow,ncol,~] = size(input);
    switch quad
        case 'L',   output = input( :, 1:floor(ncol/2), : );
        case 'R',   output = input( :, floor(ncol/2)+1:end, : );

        case 'T',   output = input( 1:floor(nrow/2), :, : );
        case 'B',   output = input( floor(nrow/2)+1:end, :, : );

        case 'TL',  output = input( 1:floor(nrow/2), 1:floor(ncol/2), : );
        case 'TR',  output = input( 1:floor(nrow/2), floor(ncol/2)+1:end, : );
        case 'BL',  output = input( floor(nrow/2)+1:end, 1:floor(ncol/2), : );
        case 'BR',  output = input( floor(nrow/2)+1:end, floor(ncol/2)+1:end, : );

        otherwise, error('Invalid subfield string %s',quad);
    end
end

output = double(output);

end %END FUNCTION subfield2


