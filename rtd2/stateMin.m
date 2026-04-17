function [] = stateMin(traceFile,dwtFile,state,minLength,outFileSel,outFileRem)
% split traces into productive (selected) and non-productive (rejected)
% events based on a minimum dwell time (minLength) in specified state

%   Copyright 2007-2015 Cornell University All Rights Reserved.
%   Updated 2016 St Jude Children's Research Hospital 


% Loading synchronized traces and associated idealizations
[dwellTimes,sampling,offsets,fretModel] = loadDWT(dwtFile);
orgTraces = loadTraces(traceFile);
idl = dwtToIdl(dwellTimes,offsets,orgTraces.nFrames,orgTraces.nTraces);

%% Creating ID list of molecules which reach selected criteria for 'accommodated' 
flagMat = zeros(orgTraces.nTraces,1);
for j = 1:length(dwellTimes)
    currentDwells = dwellTimes{j};

    for k = 1:size(currentDwells,1)
        dwellState = currentDwells(k,1);
        dwellLength = currentDwells(k,2);
        if (dwellState == state) && (dwellLength >= minLength)
            flagMat(j) = 1;
        end
    end
            
end

%%
% Selecting molecules
selTracesMat = TracesFret(orgTraces.nTraces,orgTraces.nFrames);
selTracesMat.fret = orgTraces.fret.*flagMat;
selTracesMat.donor = orgTraces.donor.*flagMat;
selTracesMat.acceptor = orgTraces.acceptor.*flagMat;

selTracesMat.fret(all(selTracesMat.fret == 0, 2), :) = [];
selTracesMat.donor(all(selTracesMat.donor == 0, 2), :) = [];
selTracesMat.acceptor(all(selTracesMat.acceptor == 0, 2), :) = [];
selTracesMat.traceMetadata = orgTraces.traceMetadata(find(flagMat));
selTracesMat.time = orgTraces.sampling * selTracesMat.time;

% Rejecting molecules
rejTracesMat = TracesFret(orgTraces.nTraces,orgTraces.nFrames);
rejTracesMat.fret = orgTraces.fret.*~flagMat;
rejTracesMat.donor = orgTraces.donor.*~flagMat;
rejTracesMat.acceptor = orgTraces.acceptor.*~flagMat;

rejTracesMat.fret(all(rejTracesMat.fret == 0, 2), :) = [];
rejTracesMat.donor(all(rejTracesMat.donor == 0, 2), :) = [];
rejTracesMat.acceptor(all(rejTracesMat.acceptor == 0, 2), :) = [];
rejTracesMat.traceMetadata = orgTraces.traceMetadata(find(~flagMat));
rejTracesMat.time = orgTraces.sampling * rejTracesMat.time;

% Generating idealizations for associated accepted/rejected traces
selIdl = idl.*flagMat;
selIdl(all(selIdl == 0, 2), :) = [];
remIdl = idl.*~flagMat;
remIdl(all(remIdl == 0, 2), :) = [];

% Saving
saveTraces(outFileSel,selTracesMat);

[selDwt,offsets] = idlToDwt(selIdl);
[path,name,~] = fileparts(outFileSel);
outDwtSel = fullfile(path, [name '.qub.dwt']);
saveDWT(outDwtSel, selDwt, offsets, fretModel, sampling);

saveTraces(outFileRem,rejTracesMat);

[remDwt,offsets] = idlToDwt(remIdl);
[path,name,~] = fileparts(outFileRem);
outDwtRem = fullfile(path, [name '.qub.dwt']);
saveDWT(outDwtRem, remDwt, offsets, fretModel, sampling);

end
