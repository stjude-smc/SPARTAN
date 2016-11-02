function data = correctTraces(data, crosstalk, scaleFluor, indexes)
% CORRECT applies crosstalk and scaling corrections to fluorescence traces.
%
%    DATA = CORRECT(DATA, CROSSTALK, SCALING) modifies the Traces object DATA
%    in place (left hand argument is optional), applying the given CROSSTALK
%    and channel SCALING corrections in that order. Use data.recalculateFret
%    to update FRET values from corrected fluorescence traces.
%      CROSSTALK(source channel number, destination channel number, trace ID)
%      SCALING(channel number, trace ID)
%
%    ... = CORRECT(..., INDEXES) only adjustes the traces listed in the
%    vector INDEXES. CROSSTALK and SCALING include values for all traces.
%
%    If the final dimension of either CROSSTALK or SCALING is unity, the
%    value will be applied to all traces.

%   Copyright 2007-2015 Cornell University All Rights Reserved.

narginchk(2,4);
nargoutchk(1,1);

nFluor = numel(data.idxFluor);


% If no traceMetadata (called from gettraces), use defaults.
% FIXME handle default values listed in fileMetadata from old versions.
if ~isfield(data.traceMetadata,'crosstalk'),
    [data.traceMetadata.crosstalk]  = deal( zeros(nFluor,nFluor) );
end

if ~isfield(data.traceMetadata,'scaleFluor'),
    [data.traceMetadata.scaleFluor] = deal( ones(nFluor,1) );
end


% If no updated crosstalk/scaling values given, use current value (no change).
if isempty(crosstalk),
    crosstalk = cat(3,data.traceMetadata.crosstalk);
end
if nargin<3 || isempty(scaleFluor),
    % Convert vectors to the same, correct, orientation.
    [data.traceMetadata.scaleFluor] = to_col(data.traceMetadata.scaleFluor);
    scaleFluor = cat(2,data.traceMetadata.scaleFluor);
end


% If a single value is supplied, apply it to all traces.
if size(crosstalk,3)==1,
    crosstalk = repmat(crosstalk, [1 1 data.nTraces]);
end
if size(scaleFluor,2)==1,
    scaleFluor = repmat(scaleFluor, [1 data.nTraces]);
end


% Verify input argument sizes match
if nargin<4,
    indexes = 1:data.nTraces;
end
if islogical(indexes), indexes=find(indexes); end

if ~(max(indexes)<=data.nTraces && size(crosstalk,3)==size(scaleFluor,2) && ...
   max(indexes)<=size(crosstalk,3)),
    error('Input argument size mismatch');
end



%% Undo previous corrections, reverting to state as acquired.
chNames = data.channelNames(data.idxFluor);  %increasing wavelength order.
% nFluor = numel(data.idxFluor);

for i=1:numel(indexes)
    
    % Undo any previous scaling.
    for ch=1:nFluor,
        data.(chNames{ch})(i,:) = data.(chNames{ch})(i,:) / ...
                                       data.traceMetadata(i).scaleFluor(ch); 
    end
    
    % Undo crosstalk subtraction in reverse order (red to blue).
    for src=nFluor:-1:1,
        for dst=nFluor:-1:1,
            if src>=dst, continue; end  %only consider forward crosstalk
            
            ch1 = chNames{src};
            ch2 = chNames{dst};
            ct = data.traceMetadata(i).crosstalk(src,dst);
            data.(ch2)(i,:) = data.(ch2)(i,:) + ct*data.(ch1)(i,:);
        end
    end

end %for each trace


%% Apply new corrections
for i=to_row(indexes),
    
    % Apply crosstalk subtraction in wavelength order (blue to red).
    for src=1:nFluor,
        for dst=1:nFluor,
            if src>=dst, continue; end  %only consider forward crosstalk
            
            ch1 = chNames{src};
            ch2 = chNames{dst};
            data.(ch2)(i,:) = data.(ch2)(i,:) - crosstalk(src,dst,i) * ...
                                                            data.(ch1)(i,:);
        end
    end
    
    % Scale fluorescence to correct for unequal brightness/sensitivity
    for ch=1:nFluor,
        data.(chNames{ch})(i,:) = data.(chNames{ch})(i,:) * scaleFluor(ch,i); 
    end
    
    % Save new correction parameters in metadata
    data.traceMetadata(i).crosstalk  = crosstalk(:,:,i);
    data.traceMetadata(i).scaleFluor = scaleFluor(:,i);

end %for each trace


end %function correctTraces


