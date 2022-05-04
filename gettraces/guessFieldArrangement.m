function output = guessFieldArrangement(input, template, channels)
% guessFieldArrangement  guess channels from frame data using heuristics
%
%   FA = guessFieldArrangement(INPUT, TEMPLATE, CHANNELS) makes a
%   reasonable initial guess as to which channels are displayed based using
%   various heuristics. FA is the field arrangement matrix whose shape
%   defines how to divide image data into channels and whose values define
%   the channel order (zeros mark a channel as unused).
%   - INPUT is the path to a movie file or a Movie object.
%   - TEMPLATE is the field arrangement matrix from current profile.
%   - CHANNELS is the struct array of channels from current profile.

narginchk(3,3);
nargoutchk(1,1);
output = [];

% Process input arguments
if nargin<1
    input = getFile('*.tif');
end
if ischar(input)
    input = Movie.load(input);
end
if isempty(input), return; end
assert( isa(input,'Movie'), 'input must be Movie object or path to movie file' );

% Index into template of the most common donor dyes (e.g., Cy3).
[~,idxGreen] = min(  abs([channels.wavelength]-532)  );

% Look for a pattern left by the Hamamatsu hard disk recorder left at
% the left (center row) of every image. If these images are tiled by the
% acquisition software (as it is in FlashGordon), we can infer the number
% and positions of the original frames from this pattern.
% The center row will start with: 0 65535 0 65535.
frame = input.readFrames(1);
low  = frame==0;
high = frame==65535;
[idxRow,idxCol] = find( low(:,1:end-3) & high(:,2:end-2) & low(:,3:end-1) & high(:,4:end) );

if numel(idxRow)>0
    % Convert marker positions into a field arrangement mask.
    r = unique(idxRow);
    c = unique(idxCol);
    output = false( numel(r), numel(c) );

    for i=1:numel(idxRow)
        % Convert pixel location of marker to field (camera) position in mosaic.
        output( idxRow(i)==r, idxCol(i)==c ) = true;
    end

    % Guess channel assignments from template
    % If the movie doesn't match the template, arbitrarily select fields
    % starting from top-left.
    % For 2-color, this will leave the top row (Cy3,Cy5).
    template = template( 1:size(output,1), 1:size(output,2) );
    template(output==false) = 0;
    output = template;

elseif input.nX >= 2*input.nY
    % If DCIMG marker not found, guess 2-color for assymetric frame size.
    output = [idxGreen idxGreen+1];
    
else
    % If all else fails, guess single channel Cy3.
    output = idxGreen;
end

end  %function

