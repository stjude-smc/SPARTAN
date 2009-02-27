function [in_a,in_b] = compareDatasets( filename_a, filename_b )
% CompareDatasets   Compare .txt traces files
% 
% [in_a,in_b] = CompareDatasets( filename_a, filename_b )
% 
% Creates two data files
% Where in_a is a list of molecular names (tags) that were found in a but
% not b, and in_b are those only in b.  Does not compare the actual data!


% Get a list of molecule names (tags)
all_a = GetTags(filename_a);
all_b = GetTags(filename_b);

in_a = setdiff(all_a,all_b); %a but not b
in_b = setdiff(all_b,all_a); %b but not a

data_a = GetData(filename_a,in_a);
data_b = GetData(filename_b,in_b);

% Save the data to disk
timeaxis = 1:1501;

added_filename = strrep(filename_a,'.txt','_added.txt');
lost_filename  = strrep(filename_a,'.txt','_lost.txt');

fid = fopen(added_filename,'w');
fprintf(fid,'%d ', timeaxis);
fprintf(fid,'\n');
fprintf(fid,'%s\n',data_a{:});
fclose(fid);

fid = fopen(lost_filename,'w');
fprintf(fid,'%d ', timeaxis);
fprintf(fid,'\n');
fprintf(fid,'%s\n',data_b{:});
fclose(fid);


end


function tags = GetTags( filename )

tags = cell(0);

fid = fopen(filename,'r');
fgetl(fid); %skip header line (time axis)


% Read in the tag names in first file
while 1,
    line = fgetl(fid);
    if ~ischar(line), break, end
    
    tag = line( 1:find(line==' ',1,'first')-1 );
    if numel(tag) < 1,  break; end
    if tag(end) == '-', tag = tag(1:end-1); end
    tags{end+1} = tag;
    
    % Skip acceptor and FRET lines (same tags)
    fgetl(fid); fgetl(fid);
    
end %for each line

fclose(fid);

end



function data = GetData( filename, tags )

data = {};

% Open input file
fid = fopen(filename,'r');
fgetl(fid); %skip header line (time axis)

% Read in the tag names in first file
while 1,
    line = fgetl(fid);
    if ~ischar(line), break, end
    
    tag = line( 1:find(line==' ',1,'first')-1 );
    if numel(tag) < 1,  break; end
    if tag(end) == '-', tag = tag(1:end-1); end
    
    i = strmatch(tag,tags);
    if numel(i)>0,
        data{end+1} = line;
    end    
    
end %for each line

fclose(fid);

end
