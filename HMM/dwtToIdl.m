function idl = dwtToIdl( dwt, traceLen, offsets )
% 

% nTraces = numel(dwt);
% 
% 
% idl = zeros(nTraces,traceLen);
% 
% 
% for dwtID=1:nTraces
%     states = dwt{dwtID}(:,1);
%     times  = double(dwt{dwtID}(:,2));
% 
%     ends = cumsum(times)+1;
%     starts = [1; ends(1:end-1)];
%     
%     for j=1:numel(states),
%         idl(dwtID,starts(j):ends(j)) = states(j);
%     end
% end


nTraces = floor(offsets(end)/traceLen)+1;

idl = zeros(nTraces,traceLen);


for dwtID=1:numel(dwt)
    traceID = floor(offsets(dwtID)/traceLen)+1;
    
    states = dwt{dwtID}(:,1);
    times  = double(dwt{dwtID}(:,2));

    ends = cumsum(times);
    starts = [1; ends(1:end-1)+1];
    
    for j=1:numel(states),
        idl(traceID,starts(j):ends(j)) = states(j);
    end
end

