function [time,occupancy] = stateOccupancy(dwtFile,traceLength)
% plot occupancy fraction over time for all states in dwtFile

%   Copyright 2007-2015 Cornell University All Rights Reserved.


[dwt,sampling,~,fretModel] = loadDWT(dwtFile);
nStates = size(fretModel,1);

occupancy = zeros(nStates,traceLength);
time = 0:sampling:(traceLength-1)*sampling;

for i=1:length(dwt)
    currentTrace = dwt{i};
    currentPosition = 1;
    allStates = currentTrace(:,1);
    for j=1:size(currentTrace,1)
        currentState = currentTrace(j,1);
        currentLength = currentTrace(j,2);
        endPosition = currentPosition + currentLength - 1;
        occupancy(currentState,currentPosition:endPosition) = occupancy(currentState,currentPosition:endPosition) + 1;
        currentPosition = endPosition+1;
    end
end

for i=1:size(occupancy,2)
    occupancy(:,i) = occupancy(:,i)./sum(occupancy(:,i));
end

[path,name,~] = fileparts(dwtFile);
outFile = fullfile(path, [name '_stateOcc.txt']);
dlmwrite(outFile,vertcat(time,occupancy));

end
