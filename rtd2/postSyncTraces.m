function [] = postSyncTraces(minFrames,preFrames,totalFrames,fileList)
% POSTSYNCTRACES Post-synchronize traces.
%
%   Loads one or more .traces or .rawtraces files and post-synchronizes the
%   traces to the first frame of non-zero FRET. The point of post-
%   synchronization is determined by SKM idealization using a simple
%   two-state model.
%
%   POSTSYNCTRACES(MINFRAMES,PREFRAMES,TOTALFRAMES) displays a dialog to
%   select trace files to be post-synchronized. Only non-zero dwells with a
%   minimum number of MINFRAMES frames are considered. A number PREFRAMES
%   of frames preceding the transition are included in the post-
%   synchronized traces. The post-synchronized traces from each original
%   file are saved (with a total length of TOTALLENGTH frames) to a file
%   with '_postSync' appended to the original file name.
%
%   POSTSYNCTRACES(MINFRAMES,PREFRAMES,TOTALFRAMES, FILELIST) uses the
%   trace files specified in the cell array FILELIST instead of prompting
%   the user.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


%%

if nargin < 4
    %ask user to select files
    fileList = getFiles('*.traces;*.rawtraces','Select trace files');
end

%parameters for skm
skmParams.maxItr = 10;
skmParams.convLL = 0.01;
skmParams.zeroEnd = 1;
skmParams.seperately = 1;
skmParams.quiet = 1;

model = QubModel(2);
model.mu    = [0 0.3];
model.sigma = [0.053 0.074];
model.p0    = [0.999 0.001];
model.rates = [0 20; 20 0];
model.fixMu(:)    = true;
model.fixSigma(:) = true;

for i = 1:numel(fileList)
    currentTraces = loadTraces(fileList{i});
    
    %idealize with two-state model and save to dwt file
    idl = skm(currentTraces.fret,currentTraces.sampling,model,skmParams);
    dwellTimes = idlToDwt(idl);
    
    [filePath,fileName] = fileparts(fileList{i});
    offsets = currentTraces.nFrames*((1:currentTraces.nTraces)-1);
    saveDWT([filePath filesep fileName '_2state.qub.dwt'],dwellTimes, ...
        offsets, model, currentTraces.sampling);
      
    % Finding number of events which satisfy criteria
    c= 0;
    for tracenum = 1:currentTraces.nTraces
        %%
        if numel(dwellTimes{tracenum}(:,1)) > 1
            state2LongDwells = find(dwellTimes{tracenum}(2:end,1) == 2 & (dwellTimes{tracenum}(2:end,2) >= minFrames)) + 1;
        else
            continue
        end

        c = c+numel(state2LongDwells);
    end

    % Pre-allocate for speed

    syncTraces = TracesFret(c,totalFrames);
    syncTraces.time = currentTraces.sampling * syncTraces.time;
    syncTraces.fileMetadata = currentTraces.fileMetadata;
    syncTraces.fret = zeros(c,totalFrames);
    syncTraces.donor = zeros(c,totalFrames);
    syncTraces.acceptor = zeros(c,totalFrames);

    c = 1;
    for j = 1:currentTraces.nTraces
        currentDwells = dwellTimes{j};

        dwellOffset = 1;
        for k = 2:size(currentDwells,1)
            dwellState = currentDwells(k,1);
            dwellLength = currentDwells(k,2);
            previousLength = currentDwells(k-1,2);
            dwellOffset = dwellOffset + previousLength;

            % select "on" dwells with specified minimum length
            if (dwellState==2) && (dwellLength >= minFrames)
                newFret = zeros(1,totalFrames);
                newDonor = zeros(1,totalFrames);
                newAcceptor = zeros(1,totalFrames);

                preLength = min(previousLength,preFrames);
                dwellLength = min(dwellLength,(totalFrames-preFrames));

                newFret((preFrames+1-preLength):(preFrames+dwellLength)) ...
                    = currentTraces.fret(j,(dwellOffset-preLength):(dwellOffset+dwellLength-1));

                newDonor((preFrames+1-preLength):(preFrames+dwellLength)) ...
                    = currentTraces.donor(j,(dwellOffset-preLength):(dwellOffset+dwellLength-1));

                newAcceptor((preFrames+1-preLength):(preFrames+dwellLength)) ...
                    = currentTraces.acceptor(j,(dwellOffset-preLength):(dwellOffset+dwellLength-1));

                % add dwell to Traces object
                syncTraces.fret(c,:) = newFret;
                syncTraces.donor(c,:) = newDonor;
                syncTraces.acceptor(c,:) = newAcceptor;
                c = c+1;
            end
        end
    end

    saveTraces([filePath filesep fileName '_postSync.traces'],syncTraces);
end

end

