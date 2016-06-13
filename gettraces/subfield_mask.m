function mask = subfield_mask(input, quad)
%SUBFIELD_MASK  Extract image subfields designated by a string
%
%   OUT = subfield(IN,STR) extracts a sub-fields from the image matrix IN
%   as described in STR, which can be empty (whole field), L (left), R (right),
%   T (top), B (bottom), TL, TR, BL, or BR. If STR is a cell array, OUT will
%   be a cell array with a subfield associated with each cell in STR.
%   
%   OUT = subfield(IN,GEO) extracts standard set of sub-fields. GEO can be:
%   1 (whole field), 2 (L/R), 3 (T/B), or 4 (quadrants).

%   Copyright 2016 Cornell University All Rights Reserved.


[nrow,ncol,~] = size(input);
mask = false(size(input));

switch upper(quad)
    case '',    mask(:)=true;  %whole field
    
    case 'L',   mask( :, 1:floor(ncol/2),     : ) = true;
    case 'R',   mask( :, floor(ncol/2)+1:end, : ) = true;
    
    case 'T',   mask( 1:floor(nrow/2),     :, : ) = true;
    case 'B',   mask( floor(nrow/2)+1:end, :, : ) = true;
    
    case 'TL',  mask( 1:floor(nrow/2),     1:floor(ncol/2),     : ) = true;
    case 'TR',  mask( 1:floor(nrow/2),     floor(ncol/2)+1:end, : ) = true;
    case 'BL',  mask( floor(nrow/2)+1:end, 1:floor(ncol/2),     : ) = true;
    case 'BR',  mask( floor(nrow/2)+1:end, floor(ncol/2)+1:end, : ) = true;
    
    otherwise, error('Invalid subfield string %s',quad);
end

assert( islogical(mask) );

end %END FUNCTION subfield_mask

