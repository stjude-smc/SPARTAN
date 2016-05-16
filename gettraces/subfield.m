function output = subfield(input, quad)
%SUBFIELD  Extract image subfields designated by a string
%
%   OUT = subfield(IN,STR) extracts a sub-fields from the image matrix IN
%   as described in STR, which can be empty (whole field), L (left), R (right),
%   T (top), B (bottom), TL, TR, BL, or BR. If STR is a cell array, OUT will
%   be a cell array with a subfield associated with each cell in STR.
%   
%   OUT = subfield(IN,GEO) extracts standard set of sub-fields. GEO can be:
%   1 (whole field), 2 (L/R), 3 (T/B), or 4 (quadrants).

%   Copyright 2016 Cornell University All Rights Reserved.


% Process input arguments
narginchk(2,2);
nargoutchk(1,1);

% Cell array of strings targeting each field to extract
if ischar(quad),
    output = subfield2(input,quad);
    
% Single string targeting a field to extract
elseif iscell(quad)
    output = cell(size(quad));
    for i=1:numel(quad)
        output{i} = subfield2(input,quad{i});
    end
    
% Integer specifing geometry. Return pre-defined set of sub-images.
elseif isnumeric(quad)
    switch quad
        case 1,  output = {input};
        case 2,  output = subfield(input, {'L','R'});
        case 3,  output = subfield(input, {'T','B'});
        case 4,  output = subfield(input, {'TL','TR','BL','BR'});
    end
else
    error('Invalid field identifier');
end

end %FUNCTION subfield



function output = subfield2(input, quad)
% Extract subfield for scalar input
[nrow,ncol,~] = size(input);

switch quad
    case '',    output = input;  %whole field
    
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

end %END FUNCTION subfield2



