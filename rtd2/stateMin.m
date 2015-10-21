function [] = stateMin(traceFile,dwtFile,state,minLength,outFileSel,outFileRem)
% split traces into productive (selected) and non-productive (rejected)
% events based on a minimum dwell time (minLength) in specified state

%   Copyright 2007-2015 Cornell University All Rights Reserved.


[dwellTimes,sampling,offsets,fretModel] = loadDWT(dwtFile);
orgTraces = loadTraces(traceFile);
idl = dwtToIdl(dwellTimes,offsets,orgTraces.nFrames,orgTraces.nTraces);

% initialize Traces object for selected traces
selTraces = TracesFret(0,orgTraces.nFrames);
selTraces.time = orgTraces.time;
selTraces.fileMetadata = orgTraces.fileMetadata;
    
remTraces = TracesFret(0,orgTraces.nFrames);
remTraces.time = orgTraces.time;
remTraces.fileMetadata = orgTraces.fileMetadata;

% initialize idl matrices
selIdl = [];
remIdl = [];

for j = 1:length(dwellTimes)
    currentDwells = dwellTimes{j};
    flag = 0;

    for k = 1:size(currentDwells,1)
        dwellState = currentDwells(k,1);
        dwellLength = currentDwells(k,2);
        if (dwellState == state) && (dwellLength >= minLength)
            flag = 1;
        end
    end
            
    if flag == 1    
        % add selected traces to Traces object and idealization
        selTraces.fret = vertcat(selTraces.fret,orgTraces.fret(j,:));
        selTraces.donor = vertcat(selTraces.donor,orgTraces.donor(j,:));
        selTraces.acceptor = vertcat(selTraces.acceptor,orgTraces.acceptor(j,:));
        selIdl = vertcat(selIdl,idl(j,:));
    else
        remTraces.fret = vertcat(remTraces.fret,orgTraces.fret(j,:));
        remTraces.donor = vertcat(remTraces.donor,orgTraces.donor(j,:));
        remTraces.acceptor = vertcat(remTraces.acceptor,orgTraces.acceptor(j,:));
        remIdl = vertcat(remIdl,idl(j,:));
    end
end

temp = TracesFret( size(selTraces.fret,1), size(selTraces.fret,2) );
temp.fret = selTraces.fret;
temp.donor = selTraces.donor;
temp.acceptor = selTraces.acceptor;
temp.fileMetadata = selTraces.fileMetadata;
temp.time = selTraces.time;

saveTraces(outFileSel,temp);

[selDwt,offsets] = idlToDwt(selIdl);
[path,name,~] = fileparts(outFileSel);
outDwtSel = fullfile(path, [name '.qub.dwt']);
saveDWT(outDwtSel, selDwt, offsets, fretModel, sampling);

temp = TracesFret( size(remTraces.fret,1), size(remTraces.fret,2) );
temp.fret = remTraces.fret;
temp.donor = remTraces.donor;
temp.acceptor = remTraces.acceptor;
temp.fileMetadata = remTraces.fileMetadata;
temp.time = remTraces.time;

saveTraces(outFileRem,temp);

[remDwt,offsets] = idlToDwt(remIdl);
[path,name,~] = fileparts(outFileRem);
outDwtRem = fullfile(path, [name '.qub.dwt']);
saveDWT(outDwtRem, remDwt, offsets, fretModel, sampling);

end
