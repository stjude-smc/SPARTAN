function idl = skmTotal( total, modelInput )
%skmTotal   Idealizes total fluorescence intensity using SKM donor blinking.
%
%   IDL = skmTotal(TOTAL) idealizes the total fluorescence intensity traces
%   in the rows of the matrix TOTAL using the segmental k-means algorithm to
%   determine when the donor is dark (total intensity is near baseline).
%   IDL is false where dark (blinking/bleached), true when it is fluorescent.
%   This method is more robust to noise traces than thresholdTotal.
%
%   ... = skmTotal(..., MODEL) uses the supplied MODEL for idealization.
%   ... = skmTotal(..., BASELINE) sets the dark state amplitude to the value.
%
%   See also: thresholdTotal, skm, TracesFret.recalculateFret.

%   FIXME: consider adding an extra parameter that defined the amplitude of
%   the "dark" state to tune the sensitivity to above-baseline noise.

%   Copyright 2015 Cornell University All Rights Reserved.

narginchk(1,2);
nargoutchk(1,1);


%% PARAMETERS
skmParams.seperately = 1;
skmParams.quiet = true;
baseline = 0.3;

% % Load a model from file, if a file name is given.
% if nargin>=2,
%     if ischar(modelInput)
%         model = qub_loadModel(modelInput);
%     elseif isscalar(modelInput) && isnumeric(modelInput),
%         baseline = modelInput;
%     else
%         error('Invalid model input');
%     end
% end

% Default model for idealization. Should this be in cascadeConstants?
% if ~exist('model','var'),
    model.class = [1 2 1]';
    model.p0    = [0.01 0.98 0.01]';
    model.mu    = [baseline 1];
    model.sigma = [0.3 0.3];
    model.rates = [0        20.0000    6.6000
                   6.6000         0         0
                   0.0001         0         0]/0.1;
    model.nStates = numel(model.class);
    model.nClass  = numel(model.mu);
% end

model.fixMu = [1 0];


%% ALGORITHM
[nTraces,nFrames] = size(total);

% Use the threshold method to get a crude estimate of total intensities
% to normalize the data to fit the model (range is 0 to 1).
idl = thresholdTotal(total);
t = sum(idl.*total,2) ./ sum(idl,2);
t(t==0) = median(t);
normTotal = bsxfun( @rdivide, total, t );

% run SKM to determine where the donor is dark.
% FIXME: change skm to return an idealization to simplify this code and
% make it faster (skm is about as fast as dwtToIdl!).
[dwt,~,~,offsets] = skm( normTotal, 100, model, skmParams );
idl = dwtToIdl(dwt, offsets, nFrames, nTraces);

% Re-normalize using the first pass idealiation and run again to refine.
t = sum(idl.*total,2) ./ sum(idl,2);
normTotal = bsxfun( @rdivide, total, t );

[dwt,~,~,offsets] = skm( normTotal, 100, model, skmParams );
idl = dwtToIdl(dwt, offsets, nFrames, nTraces);

% Remove rising/falling edges that may have low SNR and inaccurate FRET.
%idl = idl>1;
idl = imerode(idl>1,[1 1 1]);


end %function skmTotal
    
    