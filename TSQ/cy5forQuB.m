function cy5forQuB( files )
% cy5forQuB    Save Cy5 intensity data for import in QuB
%
%    
% 

const = cascadeConstants;
const.gamma=1;

% Get filename from user if not specified
if nargin < 1,
    files = getFiles;
end

nFiles = numel(files);

for i=1:nFiles,

    filename = files{i};
    
    % Load FRET data file
    data = loadTraces( filename );
    d = data.donor;
    a = data.acceptor;
    
    % Use Cy3 if channels are reversed (instead of making a Cy3forqub script).
    if mean(d(:)) > mean(a(:)), a = d; end
    d = zeros(size(a));

    % Normalize the intensity to the average.
    % Using traceStat seems to bias toward high intensities and gives a poor
    % normalization. Fitting to Gaussians sometimes works and may be more
    % accurate, but not all data has a clear Gaussian distribution.
    %stats = traceStat( d,a,f, const );
    stats = traceStat( data );
    a = a./mean([stats.t]);

%     [amp,x] = hist( a(:), 100 );
%     f0 = fit( x', amp', 'gauss2' );
%     normFactor = max( [f0.b1 f0.b2] );
% 
%     a = a./normFactor;


    % Save data to file.
    % FORMAT:
    %   All datapoints are concatinated into a M*N column vector;
    %   each datapoint is on a new line.
    %
    [p f] = fileparts(filename);
    outfilename = fullfile( p, [f '.qub_cy5.txt'] );
    
    fid=fopen(outfilename,'w');
    if fid<=0,
        error('Could not save output file (might be in use?): %s',outfilename);
    end
        
    disp( ['Saving to ' outfilename] );

    a = a';
    fprintf( fid, '%d\n', a(:) );

    fclose(fid);

end

% END FUNCTION SAVETRACESQUB

