function [dwellc] = correctedDwelltimes( dwtFilename )


% Load the dwell times
[dwells,sampling,offsets,model] = loadDWT(dwtFilename);

nTraces  = numel(dwells);
nClasses = numel(model)/2;

% Build list of dwell times in each class
dwellc    = cell(nClasses,1);

for i=1:nTraces,
    
    classes = dwells{i}(:,1);
    times   = double( dwells{i}(:,2).*sampling );
    nDwells = numel(classes);
    
    start = find(classes~=1,1,'first'); % Skip initial dwells in dark state
    curClass = classes(start);  %class last observed
    curTime  = 0;               %total time so far spent in curState
    
    for j=start:nDwells
        
        % Ignore dwells in the dark states
        if classes(j) == 1,
            curTime = curTime + times(j);
        
        % If we are still in the same state, continue summing dwell times
        elseif classes(j) == curClass,
            curTime = curTime + times(j);
        
        % If a new state is encountered, add total time in the previous
        % state and add it to the distribution.
        else
            dwellc{curClass} = [ ...
                dwellc{curClass} ; curTime ];
            
            % Reset to new internal state
            curClass = classes(j);
            curTime  = times(j);
        end
        
    end %for each dwell
        
end %for each trace

for i=1:nClasses,
    disp( numel(dwellc{i}) );
end
