function resizeTraces( traceLen )
% RESIZEQUB   Change length of a .traces file
% 
% If traceLen < actual, the traces will be cropped.
% If traceLen > actual, the last data value will be appended
% 

direct=uigetdir('','Choose directory of traces:');
if direct==0, return; end

files = dir([direct filesep '*.traces']);


constants = cascadeConstants();


% For each file in the user-selected directory
for file = files',
    
    % ---- Read traces file
    filename = [direct filesep file.name];
    [donor,accep] = LoadTraces( filename,constants );
    [nTraces,actualLen] = size(donor);
    
    % ---- Expand traces    
    if traceLen == actualLen,
        %no changes needed
        continue;
    
    elseif traceLen > actualLen,
        % Expand traces using last value
        delta = traceLen-actualLen;
        donor = [donor repmat( donor(:,end), 1, delta )];
        accep = [accep repmat( accep(:,end), 1, delta )];
    
    else
        % Truncate tracesa
        donor = donor(:,1:traceLen);
        accep = accep(:,1:traceLen);
    end
    
    disp( sprintf('Resizing %.0f to %.0f: %s',actualLen,traceLen,file.name) );
    
    
    % ---- Reconstruct traces matrix
    
    traces = zeros( nTraces*2, traceLen );
    traces(1:2:end,:) = donor;
    traces(2:2:end,:) = accep;
    
    % ---- Save modified traces to file.
    
    [p,name,e,v] = fileparts(filename);

    fid=fopen(filename,'w');

    fwrite(fid,traceLen,'int32');
    fwrite(fid,nTraces*2,'int16');  % max of 32,000 traces!

    for j=1:nTraces;
        fprintf(fid, '%s_%d-', name, j);
    end
    
    fwrite(fid,traces,'int16');
    fclose(fid);
    clear traces; 
    
end

disp('Finished.');



