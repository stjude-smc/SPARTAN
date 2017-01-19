function data = correctTraces(data, crosstalk, scaleFluor, indexes)
% correctTraces  Apply crosstalk and scaling corrections to fluorescence traces.
%
%    DATA = correctTraces(DATA, CROSSTALK, SCALING) modifies the Traces object DATA
%    in place (left hand argument is optional), applying the given CROSSTALK
%    and channel SCALING corrections in that order. Use data.recalculateFret
%    to update FRET values from corrected fluorescence traces.
%      CROSSTALK(source channel number, destination channel number, trace ID)
%      SCALING(channel number, trace ID)
%
%    ... = correctTraces(..., INDEXES) only adjustes the traces listed in the
%    vector INDEXES. CROSSTALK and SCALING include values for all traces.
%
%    If the final dimension of either CROSSTALK or SCALING is unity, the
%    value will be applied to all traces.
%
%    This function does not alter FRET traces. Use the TracesFret(4) method
%    recalculateFret to update FRET traces.
%
%    See also crosstalkcorrect, gammacorrect, scaleAcceptor, bgsub.

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

% Undo previous scaling by dividing fluorescence by the factor in metadata.
for ch=1:nFluor,
    chscale = arrayfun( @(x)x.scaleFluor(ch), data.traceMetadata(indexes) );
    data.(chNames{ch})(indexes,:) = bsxfun(@rdivide, data.(chNames{ch})(indexes,:), chscale);  %indexes
end

% Determine reverse crosstalk channel pairs.
[src,dst] = find( triu(ones(nFluor))-eye(nFluor) );
src=flip(src); dst=flip(dst);

% Undo previous crosstalk corrections by adding forward fluorescence
% by the factor in metadata.
for c=1:numel(src)
    ch1 = chNames{src(c)};
    ch2 = chNames{dst(c)};
    
    ct = arrayfun( @(x)x.crosstalk(src(c),dst(c)), data.traceMetadata(indexes) );
    data.(ch2)(indexes,:) = data.(ch2)(indexes,:) + bsxfun(@times, data.(ch1)(indexes,:), ct);
end



%% Apply new corrections

% Determine foward crosstalk channel pairs.
[src,dst] = find( triu(ones(nFluor))-eye(nFluor) );

% Subtract forward crosstalk (in wavelength order: blue to red).
for c=1:numel(src)
    ch1 = chNames{src(c)};
    ch2 = chNames{dst(c)};
    
    ct = squeeze( crosstalk(src(c),dst(c),indexes) );
    data.(ch2)(indexes,:) = data.(ch2)(indexes,:) - bsxfun(@times, ...
                                                   data.(ch1)(indexes,:), ct);
end

% Scale fluorescence to correct for unequal brightness/sensitivity
for ch=1:nFluor,
    data.(chNames{ch})(indexes,:) = bsxfun(@times, ...
                      data.(chNames{ch})(indexes,:), scaleFluor(ch,indexes)');
end

% Save new correction parameters in metadata
for i=to_row(indexes),
    data.traceMetadata(i).crosstalk  = crosstalk(:,:,i);
    data.traceMetadata(i).scaleFluor = scaleFluor(:,i);
end



end %function correctTraces


