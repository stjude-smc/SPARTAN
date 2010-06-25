function cy5forQuB( filename )
% cy5forQuB    Save Cy5 intensity data for import in QuB
%
%    
% 


% Get filename from user if not specified
if nargin < 1,
    [datafile,datapath] = uigetfile({'*.txt'},'Choose a traces.txt file:');
    if datafile == 0,
        disp('No files selected, exiting.');
        return;
    end
    
    filename = [datapath datafile];
end


% Load FRET data file
[d,a] = loadTraces( filename );

% Normalize the intensity to the average.
% Using traceStat seems to bias toward high intensities and gives a poor
% normalization. Instead I fit the intensity distribution to two Gaussians - one
% for ON and one for OFF.
% stats = traceStat( d,a,f );
% a = a./mean([stats.t]);

[amp,x] = hist( a(:), 100 );
f0 = fit( x', amp', 'gauss2' );
normFactor = max( [f0.b1 f0.b2] );

a = a./normFactor;



outfilename = strrep( filename, '.txt', '.qub_cy5.txt' );

% FORMAT:
%   All datapoints are concatinated into a M*N column vector;
%   each datapoint is on a new line.
%


a = a';

fid=fopen(outfilename,'w');
disp( ['Saving to ' outfilename] );

fprintf( fid, '%d\n', a(:) );

fclose(fid);

% END FUNCTION SAVETRACESQUB

