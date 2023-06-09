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
suffix = {'A','B','C','D'};  %add to split file names
nX = 1152;
nY = 1152;  %defaults for 2x2 binned, full-frame Hamamatsu Fusion cameras.


% Check input arguments
if nargin<1
    files = getFiles();
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
    
    % Split FOV into four quadrants by molecule position.
    % Origin is top-lelt corner. order=[A B; C D]
    A = selectCenter(data,  y <= floor(nY/2) & x <= floor(nX/2) ); % top left
    B = selectCenter(data,  y <= floor(nY/2) & x >  floor(nX/2) ); % top right
    C = selectCenter(data,  y >  floor(nY/2) & x <= floor(nX/2) ); % bottom left
    D = selectCenter(data,  y >  floor(nY/2) & x >  floor(nX/2) ); % bottom right
    
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



function output = selectCenter(data, quadrant)
% Select only the center of printed spots

STD = 1.75;  %max distance from center of printed spot

x = [data.traceMetadata(quadrant).donor_x];
y = [data.traceMetadata(quadrant).donor_y];
data = data.getSubset(quadrant);

pdx = fitdist(x','Normal');
pdy = fitdist(y','Normal');

% Define a region that is within some N standard deviations of the
% center of the density of molecules.
selected = x < pdx.mu + STD*pdx.sigma  &  x > pdx.mu - STD*pdx.sigma  ...
         & y < pdy.mu + STD*pdy.sigma  &  y > pdy.mu - STD*pdy.sigma;
output = data.getSubset(selected);

fprintf('Selected %d of %d molecules (%0.f%%)\n\n', sum(selected), ...
        numel(selected), 100*sum(selected)/numel(selected) );

% % Display of location histograms for testing.
% figure;
% subplot( numel(files), 3, (i-1)*3+1 );
% histfit(x);
% xlim([0 600]);
% subplot( numel(files), 3, (i-1)*3+2 );
% histfit(y);
% xlim([0 600]);
% subplot( numel(files), 3, (i-1)*3+3 );
% scatter( x(selected),  y(selected),  'bo' ); hold on;
% scatter( x(~selected), y(~selected), 'ro' );


end

