function forQuB( files )
% FORQUB.m converts traces files from autotrace or sorttraces into a
% format that can be imported into QuB. If multiple files are selected, data
% are combined from each file and saved to a single, merged .qub.txt file.
% To convert each trace seperately, use forQuB2.


% Get file names from user (all at once to save time)
if nargin<1,
    files = getFiles([],'Select traces files to convert');
end

% Load the traces files
fret = [];
for i=1:length(files),
    [d,a,f] = loadTraces(files{i});
    fret = [fret ; f];
end


% Create or get an output filename
[p,n] = fileparts(files{i});
outfile = [p filesep n '.qub.txt']

if length(files)>1
    outname=0;
    while outname==0,
        [outname outpath]=uiputfile(outfile,'Save new file as:');
        disp(outname);
    end
    outfile=strcat(outpath,outname);
end


% Save the data to file
data = fret';
data( data<-0.5 ) = -0.5;
data( data>1 ) = 1;

fid2=fopen(outfile,'w');
fprintf(fid2,'%f\n', data(:) );
fclose(fid2);


Ntraces = size(fret,1);
disp( sprintf('%d total FRET traces processed.',Ntraces) );


end

