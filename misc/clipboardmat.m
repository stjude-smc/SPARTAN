function output = clipboardmat(varargin)
%clipboardmat  Copy matrix data to clipboard
%
%   CLIPBOARDMAT(MAT) copies the matrix MAT in text form into the clipboard.
%   The clipboard('copy',...) function only supports text.
%
%   CLIPBOARDMAT(hObject,eventData,MAT) can be used as a callback function:
%     uimenu('Label','Copy', 'Callback',{@clipboardmat,output});

%   Copyright 2016 Cornell University All Rights Reserved.


nargoutchk(0,1);

switch nargin
    case 1
        input = varargin{1};
    case 3  %used as a callback with handle inputs first
        input = varargin{3};
    otherwise
        error('Invalid number of arguments');
end

output = sprintf(['%f' repmat(' %f',1,size(input,2)-1) '\n'], input');
clipboard('copy',output);

end