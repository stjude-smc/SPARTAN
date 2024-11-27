function selectPrintedSpots( files )
% Select molecules with positions that are within four quadrants of the
% field of view multiplexed surface patterning experiments. To reduce the
% contribution of non-specifically bound molecules, only the center of each
% spot is selected using a 2D Gaussian fit of molecule locations.
% 
% Files will be saved with the following extensions:
%  - A (top left)
%  - B (top right)
%  - C (bottom left)
%  - D (bottom right)

% Parameters
STD = 1.9;  %max distance from center of printed spot, number of standard deviations
suffix = {'A','B','C','D'};  %add to split file names
nX = 1152;
nY = 1152;  %defaults for 2x2 binned, full-frame Hamamatsu Fusion cameras.


% Check input arguments
if nargin<1
    files = getFiles('.rawtraces');
elseif ischar(files)
    files = {files};
elseif ~iscell(files)
    error('Input must be a string or cell array of strings');
end



for i=1:numel(files)
    
    % Load trace data and extract molecule locations
    data = loadTraces(files{i});
    [p,f,e] = fileparts(files{i});
    x = to_row( [data.traceMetadata.donor_x] );
    y = to_row( [data.traceMetadata.donor_y] );
    
    if isfield(data.fileMetadata,'nX')
        nX = data.fileMetadata.nX;
        nY = data.fileMetadata.nY;
    end
    
    figure; hold on;
    title([f e],'Interpreter','none');
    %legend( {'Rejected','A','B','C','D'} );
    axis([0 nX 0 nY]);
    
    % Split FOV into four quadrants by molecule position.
    % Origin is top-lelt corner. order=[A B; C D]
    A = selectCenter(data,  y <= floor(nY/2) & x <= floor(nX/2), STD, 'ro' ); % top left
    B = selectCenter(data,  y <= floor(nY/2) & x >  floor(nX/2), STD, 'bo' ); % top right
    C = selectCenter(data,  y >  floor(nY/2) & x <= floor(nX/2), STD, 'go' ); % bottom left
    D = selectCenter(data,  y >  floor(nY/2) & x >  floor(nX/2), STD, 'mo' ); % bottom right
    
    % Save the each subset associated with a quadrant
    outname = fullfile( p, [f '_' suffix{1} e] );
    saveTraces( outname, A );
    
    outname = fullfile( p, [f '_' suffix{2} e] );
    saveTraces( outname, B );
    
    outname = fullfile( p, [f '_' suffix{3} e] );
    saveTraces( outname, C );
    
    outname = fullfile( p, [f '_' suffix{4} e] );
    saveTraces( outname, D );

end %for each file



end %function



function output = selectCenter(data, quadrant, STD, color)
% Select only the center of printed spots

x = [data.traceMetadata(quadrant).donor_x];
y = [data.traceMetadata(quadrant).donor_y];
data = data.getSubset(quadrant);

pdx = fitdist(x','Normal');
pdy = fitdist(y','Normal');


% Define a region that is within some N standard deviations of the
% center of the density of molecules.
selected = ((x - pdx.mu)/pdx.sigma).^2 + ((y - pdy.mu)/pdy.sigma).^2 < STD.^2;

output = data.getSubset(selected);

fprintf('Selected %d of %d molecules (%0.f%%)\n\n', sum(selected), ...
        numel(selected), 100*sum(selected)/numel(selected) );

scatter( x(selected),  y(selected),  color );
scatter( x(~selected), y(~selected), 'ko' );

end

