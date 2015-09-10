function [] = simplePostSync(preFrames,totalFrames,fretThreshold,fileList)
% like postSyncTraces, but uses non-iterative, threshold-based
% post-synchronization

%   Copyright 2007-2015 Cornell University All Rights Reserved.


if nargin < 4
    %ask user to select files
    fileList = getFiles('*.traces;*.rawtraces','Select trace files');
end

nFiles = length(fileList);

for i = 1:nFiles
    currentTraces = loadTraces(fileList{i});
    [filePath,fileName,fileExt] = fileparts(fileList{i});
      
    % initialize Traces object
    syncTraces = TracesFret(0,totalFrames);
    syncTraces.time = currentTraces.sampling * syncTraces.time;
    syncTraces.fileMetadata = currentTraces.fileMetadata;
    
    for j = 1:currentTraces.nTraces
        currentFret = currentTraces.fret(j,:);
        currentDonor = currentTraces.donor(j,:);
        currentAcceptor = currentTraces.acceptor(j,:);
        
        newFret = zeros(1,totalFrames);
        newDonor = zeros(1,totalFrames);
        newAcceptor = zeros(1,totalFrames);
        
        previousFrames = find(currentFret >= fretThreshold,1,'first') - 1;
        preLength = min(previousFrames,preFrames);
        remainingLength = min(size(currentFret,2)-previousFrames,totalFrames-preFrames);
        
        newFret((preFrames+1-preLength):(preFrames-preLength+remainingLength)) ...
            = currentTraces.fret(j,(previousFrames+1-preLength):(previousFrames-preLength+remainingLength));
                
        newDonor((preFrames+1-preLength):(preFrames-preLength+remainingLength)) ...
            = currentTraces.donor(j,(previousFrames+1-preLength):(previousFrames-preLength+remainingLength));

        newAcceptor((preFrames+1-preLength):(preFrames-preLength+remainingLength)) ...
            = currentTraces.acceptor(j,(previousFrames+1-preLength):(previousFrames-preLength+remainingLength));
        

        % add dwell to Traces object
        syncTraces.fret = vertcat(syncTraces.fret,newFret);
        syncTraces.donor = vertcat(syncTraces.donor,newDonor);
        syncTraces.acceptor = vertcat(syncTraces.acceptor,newAcceptor);
    end
    temp = TracesFret( size(syncTraces.fret,1), size(syncTraces.fret,2) );
    temp.fret = syncTraces.fret;
    temp.donor = syncTraces.donor;
    temp.acceptor = syncTraces.acceptor;
    temp.fileMetadata = syncTraces.fileMetadata;
    temp.time = syncTraces.time;
    
    saveTraces([filePath filesep fileName '_simplePostSync.traces'],temp);
%    outFileList{i} = [filePath filesep fileName '_simplePostSync.traces'];
%    plotTitleList{i} = fileName;
end

%plotOptions.contour_length = totalFrames;
%makeplots(outFileList,plotTitleList);%,plotOptions);

end

