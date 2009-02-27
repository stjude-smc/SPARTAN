function rebintraces(factor, files)
%REBINTRACES  Creates traces with (simulated) reduced time resolution
% 
%   REBINTRACES(FILES)
%   Shrinks the temporal resolution of trace data by FACTOR (0.1 means half
%   the original resolution)
%   If no FILES are specified, user will be asked to select them.
%   

% TODO: add varargin options for stuff

% INITIALIZE & PROCESS FUNCTION PARAMETERS

assert( factor>=1, 'Cannot expand trace resolution' );


%% INITIALIZE & PROCESS FUNCTION PARAMETERS

% If no files specified, prompt user for them.
if ~exist('files','var'),
    files = cell(0,1);

    disp('Select traces files, hit cancel when finished');
    while 1,
        [datafile,datapath] = uigetfile({'*.txt'},'Choose a traces file:');
        if datafile==0, break; end  %user hit "cancel"

        files{end+1} = [datapath filesep datafile];
    end
end

nFiles = numel(files);

if nFiles == 0,
    disp('No files specified, exiting.');
    return;
end


%% 
for i=1:nFiles,
    
    % Load traces file
    [donor,acceptor,fret,ids] = LoadTraces( files{i} );
    
    % Shrink it
    newSize = ceil(size(donor,2)/factor);
    nTraces = size(donor,1);
    disp( sprintf('Resizing to %d frames',newSize) );
    
    donor2 = zeros(nTraces,newSize);
    acceptor2 = zeros(nTraces,newSize);
    
    for j=1:nTraces,
       donor2(j,:)    = TimeScale( donor(j,:),    factor );
       acceptor2(j,:) = TimeScale( acceptor(j,:), factor );  
    end
    times2 = times(1:newSize); %*factor;
    
    fret = acceptor2./(donor2+acceptor2);
    
    
    % Write file back to disk
    filename = strrep(files{i}, '.txt', '_shrunk.txt');
    
    fid=fopen(filename,'w');

    % Write fluorescence data: {name} {datapoints...}
    % 3 lines per molecule: donor, acceptor, fret
    fprintf(fid, '%f ', times2);
    fprintf(fid,'\n');
    
    for j=1:nTraces
        
        fprintf(fid,'%s ', ids{j});
        fprintf(fid,'%g ', donor2(j,:));
        fprintf(fid,'\n');
        
        fprintf(fid,'%s ', ids{j});
        fprintf(fid,'%g ', acceptor2(j,:));
        fprintf(fid,'\n');

        fprintf(fid,'%s ', ids{j});
        fprintf(fid,'%g ', fret(j,:));
        fprintf(fid,'\n');
        
    end % for each molecule
    
    fclose(fid);
end

% end function


function output = TimeScale( trace, factor )

assert( factor >= 1 );

nFrames2 = ceil(numel(trace)/factor);
output = zeros( 1,nFrames2 );

j = 1;
itr2 = 1;
while floor(j)+1 <= numel(trace),
    
    itr1 = floor(j);
    
    r1 = 1 - (j-floor(j));
    r2 = factor-r1;
    
    output(itr2) = (r1*trace(itr1) + r2*trace(itr1+1))/factor;
    
    j = j+factor;
    itr2 = itr2+1;
    
end










