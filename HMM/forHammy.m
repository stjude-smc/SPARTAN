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

disp( 'WARNING: this program will create a file for each trace!' );
disp( 'Best to include only a few of your best traces to define the model' );


% Get filename from user
[name,filepath]=uigetfile('*.txt','Choose a fret file:');
if name==0,  return;  end
filename = strcat( filepath, filesep, name );


% Load the traces files
[donor,acceptor,fret,ids] = loadTraces( filename );
[Ntraces,Nframes] = size(donor);

% Calculate lifetimes
lifetimes = CalcLifetime( donor+acceptor )-4;


% Save data to file
basename = strrep(filename,'.txt','_HaMMy');
basename = strrep(basename,'.','p');

for i=1:Ntraces, % for each trace

    % Create output file for this trace
    outfile = sprintf( '%s%04d.dat', basename, i );
    
    % Prune data to remove donor-dark regions, which confuse HaMMy
    window = fret(i,:) ~= 0;
    
    % Gather output data in correct format
    % (truncate to before photobleaching)
    len = sum(window);
    data = [1:len ; donor(i,window) ; acceptor(i,window)];
    
    % Write file to disk as a vector (columnwise through <data>)
    dlmwrite( outfile, data(:)', ' ' );

end % for each trace


disp( sprintf('Sucessfully saved %d FRET traces', Ntraces) );

end % function forHammy

