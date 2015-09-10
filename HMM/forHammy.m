function forHammy()
% FORHAMMY  Converts fluorescence data into HaMMy format
%
%   Prompts user for location of traces file (output from autotrace)
%   to convert.  Output extension is .dat.  The original is left intact.
%   
%   REQUIRES: LoadTraces.m
%
% http://bio.physics.uiuc.edu/HaMMy.html
% FORMAT: <time> <donor intensity> <acceptor intensity> ...

%   Copyright 2007-2015 Cornell University All Rights Reserved.


disp( 'WARNING: this program will create a file for each trace!' );
disp( 'Best to include only a few of your best traces to define the model' );


% Get filename from user
[name,filepath]=uigetfile('*.traces','Choose a fret file:');
if name==0,  return;  end
filename = strcat( filepath, filesep, name );


% Load the traces files
data = loadTraces( filename );
[Ntraces,Nframes] = size(data.donor);

% Calculate lifetimes
lifetimes = calcLifetime( data.donor+data.acceptor )-4;


% Save data to file
[p,f] = fileparts(filename);
basename = fullfile( p, [f '_HaMMy'] );
basename = strrep(basename,'.','p');

for i=1:Ntraces, % for each trace

    % Create output file for this trace
    outfile = sprintf( '%s%04d.dat', basename, i );
    
    % Prune data to remove donor-dark regions, which confuse HaMMy
    window = data.fret(i,:) ~= 0;
    
    % Gather output data in correct format
    % (truncate to before photobleaching)
    len = sum(window);
    output = [1:len ; data.donor(i,window) ; data.acceptor(i,window)];
    
    % Write file to disk as a vector (columnwise through <data>)
    dlmwrite( outfile, output(:)', ' ' );

end % for each trace


disp( sprintf('Sucessfully saved %d FRET traces', Ntraces) );

end % function forHammy

