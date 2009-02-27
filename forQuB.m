function forQuB
% FORQUB.m converts traces files from autotrace or sorttraces into a
% format that can be imported into QuB.


% Get file names from user (all at once to save time)
files=cell(0,1);
filepath = '';

while 1
    [name,filepath]=uigetfile('*.txt','Choose a fret file:',filepath);
    if name==0,  break;  end
    
    files{end+1} = strcat(filepath,name);
    
    disp('Another file? If not, press cancel.');
end

if isempty(files)
    disp('No files selected.');
    return;
end


% Load the traces files
fret = [];
for i=1:length(files),
    [d,a,f] = loadTraces(files{i});
    fret = [fret ; f];
end


% Create or get an output filename
if length(files)>1
    outname=0;
    while outname==0,
        [outname outpath]=uiputfile('*.txt','Save new file as:');
        disp(outname);
    end
    outfile=strcat(outpath,outname);
else
    outfile=strrep(files{1},'.txt','.qub.txt');
end


% Save the data to file
data = fret';

fid2=fopen(outfile,'w');
fprintf(fid2,'%f\n', data(:) );
fclose(fid2);


Ntraces = size(fret,1);
disp( sprintf('%d total FRET traces processed.',Ntraces) );


end

