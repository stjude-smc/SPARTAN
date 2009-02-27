% Catalog of FRET values and corrosponding distances from LeuT data
%



%% Load matrix of descriptors for each dataset
clear
fid = fopen('catalog.txt');
C = textscan(fid, '%q\t%q\t%q\t%q\t%*q');
nFiles = length(C{1});

% Construct a cell array for the results:
% pos1  pos2    mut.    cond.   fname   data
data = [C{:} cell(nFiles,1) cell(nFiles,1) ];


%% Get filenames from user to fill catalog

skipExisting = 1;

% List of files to analyze
prefix = '/home/dsterry/cornell/data/Daniel/Youngfang';
prev = prefix;
filter = {'*.txt','Text files (*.txt)';'*auto.txt','Autotrace Files (*auto.txt)'};
nFiles = size(data,1);

for i=1:nFiles,
    
    disp(i);
    
    % If file has already been found and tagged, skip it
    filename = data{i,5};
    if ~isempty(filename) && exist(filename,'file') && skipExisting,
        continue;
    end
    
    % Ask user for file location
    prompt = ['Select: ' data{i,1} '/' data{i,2} ' ' data{i,3} ' ' data{i,4} ':'];
    
    [filename,pathname] = uigetfile( filter, prompt, prev );
    if filename==0, continue; end
    
    fname = [pathname filename];
    data{i,5} = fname;
    
    prev = fileparts(fname);
    prev = prev( 1:find(prev=='/',1,'last') );

end

fileCatalog = data;
save( 'catalog.mat', 'fileCatalog' );




%% Fill catalog with data...

titles = cell(nFiles,1);

for i=1:nFiles,
    titles{i} = [data{i,1} '/' data{i,2} ' ' data{i,3} ' ' data{i,4} ':'];
end

params = fitPeaks( {data{:,5}}, titles );

%%

fnames = {data{:,5}};
idx = find( ~cellfun( @isempty, fnames ) );
data(idx,6) = params';


%% Extract mean values (in order) and export

for i=1:nFiles,
    model = data{i,6};
    if ~isempty(model)
        mu = sort( model([1,3]) );
        disp( mu );
    else
        disp( sprintf('    -         -') );
    end
end





















