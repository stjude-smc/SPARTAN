function forvbFRET()
% FORHAMMY  Converts fluorescence data into HaMMy format
%
%   Prompts user for location of traces file (output from autotrace)
%   to convert.  Output extension is .dat.  The original is left intact.
%   
% http://vbfret.sourceforge.net/
% FORMAT: <donor intensity> <acceptor intensity> ...
% one column per trace, two columns per molecule.
% NOTE: donor dark regions are removed before passing it on.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Get filename from user
[name,filepath]=uigetfile('*.traces','Choose a fret file:');
if name==0,  return;  end
filename = strcat( filepath, filesep, name );


% Load the traces files
data = loadTraces( filename );
[Ntraces,Nframes] = size(data.donor);

% Calculate lifetimes
% lifetimes = CalcLifetime( donor+acceptor )-2;

% Transform data into output format.
output = zeros(Nframes,Ntraces);
for i=1:Ntraces,
    % Prune data to remove donor-dark regions, which confuse HaMMy
    window = data.fret(i,:) ~= 0;
    data.donor(i,~window) = 0;
    data.acceptor(i,~window) = 0;
    
    % 
    output(:,(2*(i-1))+1) = data.donor(i,:)';
    output(:,(2*i)    ) = data.acceptor(i,:)';
end

% Save data to file
outfile = strrep(filename,'.traces','_vbFRET.dat')

% Write header (trace names).
fid = fopen(outfile,'w');
ids = num2cell( ceil(0.5:0.5:Ntraces) );
[~,baseName] = fileparts(filename);
baseName = strrep( baseName, ' ','_' );

fprintf( fid, [baseName '%04d\t' baseName '%04d\t'], ids{:} );
fprintf( fid, '\n' );

% Write file to disk as a vector (columnwise through <data>)
dlmwrite( outfile, output, '-append', 'delimiter','\t' );


end % function forvbFRET

