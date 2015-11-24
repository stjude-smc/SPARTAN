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

%   Copyright 2015 Cornell University All Rights Reserved.

%   FIXME: consider adding an extra parameter that defined the amplitude of
%   the "dark" state to tune the sensitivity to above-baseline noise.

narginchk(1,2);
nargoutchk(1,1);


%% PARAMETERS
skmParams.seperately = 1;
skmParams.quiet = true;
baseline = 0.2;

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
% The second state is a partially quenched blinking state.
% if ~exist('model','var'),
    model.class = [1 2 3]';
    model.p0    = [0.01 0.01 0.98]';
    model.mu    = [0 baseline 1];
    model.sigma = [0.15 0.15 0.3];
    model.rates = [0        0       0
                   0.1      0       2    %bleaching, ressurection rates.
                   0.1      1       0];  %blinking rate (s-1)
    model.nStates = numel(model.class);
    model.nClass  = numel(model.mu);
% end

% The blinking value is fixed so it doesn't creep up to the "on" state value.
model.fixMu = [0 1 0];


%% ALGORITHM
[nTraces,nFrames] = size(total);

% Use the threshold method to get a crude estimate of total intensities
% to normalize the data to fit the model (range is 0 to 1).
idl = thresholdTotal(total);
t = sum(idl.*total,2) ./ sum(idl,2);
t(t==0) = median(t);
normTotal = bsxfun( @rdivide, total, t );

% Run SKM to determine where the donor is dark.
[dwt,~,~,offsets] = skm( normTotal, 100, model, skmParams );
idl = dwtToIdl(dwt, offsets, nFrames, nTraces)==3;

% Re-normalize using the first pass idealiation and run again to refine.
t = sum(idl.*total,2) ./ sum(idl,2);
normTotal = bsxfun( @rdivide, total, t );

[dwt,~,~,offsets] = skm( normTotal, 100, model, skmParams );
idl = dwtToIdl(dwt, offsets, nFrames, nTraces)==3;

% Remove rising/falling edges that may have low SNR and inaccurate FRET.
idl = imerode(idl, [1 1 1]);



end %function skmTotal
    
    