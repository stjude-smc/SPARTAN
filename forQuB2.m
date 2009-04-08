function forQuB2( files )
% FORQUB2  Converts auto.txt or .traces files into .qub.txt files 
%
%   forQuB2( FILENAMES )
%   Loads each file specified and saves the FRT data to a .qub.txt file.


% Get file names from user (all at once to save time)
if nargin<1,
    files = getFiles();
end

for i=1:length(files),
    [d,a,fret] = loadTraces( files{i} );
    fret = fret';

    % Create or get an output filename
    outfile=strrep(files{i},'.traces','.txt');
    outfile=strrep(outfile,'.txt','.qub.txt');
    
    % Save the data to file
    fid2=fopen(outfile,'w');
    fprintf(fid2,'%f\n', fret(:) );
    fclose(fid2);
    
    Ntraces = size(fret,2);
    disp( sprintf('%d total FRET traces processed.',Ntraces) );
end

