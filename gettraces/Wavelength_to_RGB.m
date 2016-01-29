function RGB = Wavelength_to_RGB(wl)
% Wavelength_to_RGB: converts wavelength to RGB colors (0 to 1).
%
%     RGB = Wavelength_to_RGB(wl)
%
% wl is the wavelength of light (380 to 780 nm), RGB is a row vector of the
% RGB channel intensities (0 to 1). For example, 510 = green = [0 1 0]. If
% wl is an Nx1 vector of wavelengths, RGB is a Nx3 matrix with each
% wavelength represented in a separate row.

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% This script inspired by the method described:
% http://scslin.blogspot.com/2011/12/matlab-calculate-rbg-of-visible.html
% 


% If no parameters given, display the full spectrum
if nargin<1,
    spectrum();
    return;
end

assert( any(size(wl)==1), 'List of wavelengths must be a vector' );

RGB = zeros( numel(wl),3 ); %black outside visible range

for i=1:numel(wl),
    WL = wl(i);
    
    % Calculate RGB color values across the visible spectrum.
    B = minmax( (510-WL)/(510-490) );

    if WL<=510
        R = minmax( (440-WL)/(440-380) );
        G = minmax( (WL-440)/(490-440) );

    elseif WL>=510
        R = minmax( (WL-510)/(580-510) );
        G = minmax( (645-WL)/(645-580) );
    end

    % Fade to black at the limits of human vision (<420, >700).
    SSS = minmax( (780-WL)/(780-680) )*  minmax( (WL-380)/(420-380) );

    RGB(i,:) = SSS.*[R G B];
    
end %for each color



end % function


function val = minmax( val )
% Fix the range of a number to be [0,1]

val = max( 0, min(1,val) );

end


function spectrum()
% Display the spectrum of colors by wavelength.

RGB = Wavelength_to_RGB( 380:780 );
CV = zeros(1,380+size(RGB,1),3);
CV(1,380:780,:) = RGB;

h = figure;  cax = axes;
imagesc(CV, 'Parent',cax)
set(cax,'xlim',[380 780],'yticklabel','');
xlabel(cax,'Wavelength (nm)');
set(h,'Position',[465 564 560 94]);

end

