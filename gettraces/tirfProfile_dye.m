function [conversion,efficiency,uniformity] = tirfProfile_dye(files)
% tirfProfile_dye  Display and analysis of TIRF illumination profile
%
%  [conversion,efficiency,uniformity] = tirfProfile_dye(files);
%  files: cell array of paths to movie files.
%  conversion: multiply laser power to get kW/cm^2 illumination intensity.
%  efficiency: fraction of illumination power falling within field of view.
%  uniformity: RMS deviation of intensity from average (lower is better).
%
% This function designed for 20x magnification where approximately half of
% the pixels are illuminated.

% TO DO: add 1D profile overlay X and Y.
% TO DO: consider showing 1D profiles in a separate figure so they're easier to see.
%   or use custom panel placement like makeplots.


%% PARAMETERS

% Parameters for microscope with Fusion cameras
pxSize = (2*6.5)/20;     %Pixel size in input image (20x, Fusion cameras).
fieldX = (2304*6.5)/60;  %width of 60x FOV (microns)
fieldY = (2304*6.5)/60;  %height of 60x FOV (microns)
stepX  = 450;            %new LANE step size (microns)
stepY  = 450;            %new FIELD step size (microns)
cplot_scale = 6;         %threshold for scaling contour plot

% Uncomment these lines for narrowed field of view (no cyl. lens telescope).
% exposure = 4;  %ms, camera exposure time, which determines field size.
% fieldY = fieldY * exposure/11.22;  %11.22 is the min. full-field integration time.
% stepY = 150;
% cplot_scale = 18;

colorOrder = 'bgrk';


%% 2D contour plots and associated statistics
if nargin<1
    files = getFiles('*.tif*');
end
ax = zeros( numel(files)+1, 1 );
[conversion,efficiency,uniformity] = deal(  zeros(numel(files),1) );
if isempty(files), return; end

figure;
names = trimtitles(files);
fields = cell( numel(files), 1 );

for f=1:numel(files)
    m = Movie.load( files{f} );

    % Average together all frames
    fov = zeros(m.nY, m.nX);
    for i=1:min(m.nFrames,100)
        fov = fov + double( m.readFrames(i) );
    end

    % Subtract baseline and normalize so FOV sums to 1.
    sortpx = sort( fov(:) );
    fov = fov - mean( sortpx(1:floor(numel(sortpx)/20)) );    
    fov = fov ./ sum(fov(:));
    fields{f} = fov;
    
    % Use center of mass of illuminated area to align field of view.
    % FIXME: this really should just be the true center to match 60x?
    measurements = regionprops( true(size(fov)), fov, 'WeightedCentroid' );
    com = measurements.WeightedCentroid * pxSize;
    
    % Estimate the fraction of intensity that falls within 60x FOV.
    left   = floor( (com(1)-fieldX/2) /pxSize );
    right  = floor( (com(1)+fieldX/2) /pxSize );
    top    = floor( (com(2)-fieldY/2) /pxSize );
    bottom = floor( (com(2)+fieldY/2) /pxSize );
    inside = fov( top:bottom, left:right  );
    
    %left    = floor( (com(1)-stepX/2) /pxSize );
    %right   = floor( (com(1)+stepX/2) /pxSize );
    %top     = floor( (com(2)-stepY/2) /pxSize );
    %bottom  = floor( (com(2)+stepY/2) /pxSize );
    %outside = fov( top:bottom, left:right  );
    %disp( sum(outside(:)) / sum(fov(:)) );  %sanity check; should be ~0.
    
    efficiency(f) = sum(inside(:)) / sum(fov(:));
    conversion(f) = efficiency(f) / (fieldX*fieldY*(1e-4^2)) / 1000;
    
    % Estimate RMS uniformity relative to mean intensity w/i 60x FOV.
    mp = mean(inside(:));
    uniformity(f) = sqrt(  sum((inside(:)-mp).^2) /numel(inside)  )  /mp;
    
    % Display beam profile to scale
    ax(f) = subplot(1,numel(files)+1,f);
    imshow( fov, [0 cplot_scale/numel(fov)], ...
            'XData',[0 pxSize*size(fov,2)], 'YData',[0 pxSize*size(fov,1)] );
    colormap('jet');
    title( names{f} );
    hold on;
    %if f==numel(files), colorbar; end
    
    plot( com(1), com(2), 'k*' ); %center of mass position
    rectangle( 'Position', [ com(1)-fieldX/2 com(2)-fieldY/2 fieldX fieldY] );  %60x field of view
    rectangle( 'Position', [ com(1)-stepX/2 com(2)-stepY/2 stepX stepY] );  %mask area (move size)
    drawnow;
    
    % Display beam profile statistics
    xlabel( {sprintf('Conversion: %.3f',conversion(f)), ...
             sprintf('Efficiency: %.0f%%',efficiency(f)*100), ...
             sprintf('Uniformity: %.1f%%',uniformity(f)*100) } );
    
    % Save contour data for making figures.
%     [p,fname] = fileparts( files{f} );
%     save(  fullfile(p,[fname '_cplot.txt']), '-ASCII', 'fov' );

    % 1D profiles
    ax(end) = subplot( 1, numel(files)+1, numel(files)+1 );
    profile = sum(fov,1);
    plot( pxSize*(0:numel(profile)-1), profile, colorOrder(f) );
    
    xlabel('X Position (um)');
    hold on;
    pbaspect([1 1 1])
end

legend( ax(end), names );
% linkaxes(ax,'xy');


end  %function

