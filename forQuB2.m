function forQuB2
% FORQUB.m converts traces files from autotrace or sorttraces into a
% format that can be imported into QuB.


% Get file names from user (all at once to save time)
files = getFiles();

for i=1:length(files),
    [d,a,fret] = loadTraces(files{i});
    fret = fret';

    % Create or get an output filename
    outfile=strrep(files{i},'.txt','.qub.txt');
    
    % Save the data to file
    fid2=fopen(outfile,'w');
    fprintf(fid2,'%f\n', fret(:) );
    fclose(fid2);
    
    Ntraces = size(fret,2);
    disp( sprintf('%d total FRET traces processed.',Ntraces) );
end

