function forvbFRET()
% FORHAMMY  Converts fluorescence data into HaMMy format
%
%   Prompts user for location of traces file (output from autotrace)
%   to convert.  Output extension is .dat.  The original is left intact.
%   
%   REQUIRES: LoadTraces.m

% PROBLEM: need to remove undefined FRET regions

% http://vbfret.sourceforge.net/

% FORMAT: <donor intensity> <acceptor intensity> ...
% one column per trace, two columns per molecule.


% Get filename from user
[name,filepath]=uigetfile('*.txt','Choose a fret file:');
if name==0,  return;  end
filename = strcat( filepath, filesep, name );


% Load the traces files
[donor,acceptor] = loadTraces( filename );
[Ntraces,Nframes] = size(donor);

% Calculate lifetimes
% lifetimes = CalcLifetime( donor+acceptor )-2;

% Transform data into output format.
output = zeros(Nframes,Ntraces);
for i=1:Ntraces,
    % Prune data to remove donor-dark regions, which confuse HaMMy
    window = fret(i,:) ~= 0;
    donor(i,:) == 0 ) = 0;
    
    % Gather output data in correct format
    % (truncate to before photobleaching)
%     len = sum(window);
%     data = [1:len ; donor(i,window) ; acceptor(i,window)];
    
    output(:,(2*(i-1))+1) = donor(i,:)';
    output(:,(2*i)    ) = acceptor(i,:)';
end

% Save data to file
% Write file to disk as a vector (columnwise through <data>)
outfile = strrep(filename,'.txt','_vbFRET.dat');
dlmwrite( outfile, output, ' ' );


end % function forvbFRET

