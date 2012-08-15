function resizeTraces( traceLen )
% resizeTraces   Change length of a .traces files
%
%   resizeTraces( TRACE_LEN )
%   loads all .traces files in the directory specified by the user,
%   resizes them to the TRACE_LEN, saves them back to their original
%   filenames. This ensures that all files will be the same length.
%
%   If traceLen < actual, the traces will be cropped.
%   If traceLen > actual, the last data value will be appended
% 

if nargin<1,
    f = inputdlg('How many frames to keep?');
    if isempty(f), return; end
    traceLen = str2double(f);
end

direct=uigetdir('','Choose directory of traces:');
if direct==0, return; end

files = dir([direct filesep '*.traces']);


constants = cascadeConstants();


% For each file in the user-selected directory
for file = files',
    
    % ---- Read traces file
    filename = [direct filesep file.name];
    [donor,accep,f,ids,time] = LoadTraces( filename );
    [nTraces,actualLen] = size(donor);
    
    % ---- Reconstruct time axis (assuming linear, contiguous acquisition)
    dt = time(2)-time(1);
    time = cumsum( [time(1) repmat(dt,1,traceLen-1)] );
    
    % ---- Undo crosstalk correction (otherwise it will be done twice).
    accep = accep + constants.crosstalk*donor;
    
    % ---- Modify traces, if they do not match the target trace length.
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
    
    % ---- Save the results.
    saveTraces( filename,'traces',donor,accep,ids,time );
    
    
end

disp('Finished.');



