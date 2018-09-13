function idl = skmTotal( total, modelInput )
%skmTotal   Idealize total fluorescence intensity using SKM
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

%   Copyright 2015-2016 Cornell University All Rights Reserved.

%   FIXME: consider adding an extra parameter that defined the amplitude of
%   the "dark" state to tune the sensitivity to above-baseline noise.

narginchk(1,2);
nargoutchk(1,1);


%% PARAMETERS
skmParams.seperately = 1;
skmParams.quiet = true;
skmParams.maxItr = 20; %typically 3-5, very rarely more than 10.

% Define a default model: 1) dark, 2) blinking/quenched, 3) ON.
model = QubModel(3);
model.p0    = [0.01 0.01 0.98]';
model.mu    = [0 0.2 1];   %second value defines sensitivity to partially quenched states.
model.sigma = [0.15 0.15 0.3];
model.rates = [0        0       0
               0.1      0       2    %bleaching, ressurection rates.
               0.1      1       0];  %bleaching, blinking rate (s-1)

% Use a user-defined model instead
if nargin>=2
    if isa(modelInput,'QubModel')
        model = modelInput;
    elseif ischar(modelInput)
        model = QubModel(modelInput);
    elseif isnumeric(modelInput) && numel(modelInput)==1
        model.mu(2) = modelInput;
    else
        error('Invalid model input')
    end
end

% The blinking value is fixed so it doesn't creep up to the "on" state value.
model.fixMu = [0 1 0];


%% ALGORITHM

% Use the threshold method to get a crude estimate of total intensities
% to normalize the data to fit the model (range is 0 to 1).
idl = thresholdTotal(total);
t = sum(idl.*total,2) ./ sum(idl,2);
t(t==0) = median(t);
normTotal = bsxfun( @rdivide, total, t );

% Run SKM to determine where the donor is dark.
idl = skm( normTotal, 100, model, skmParams )==3;

% Re-normalize using the first pass idealiation and run again to refine.
t = sum(idl.*total,2) ./ sum(idl,2);
normTotal = bsxfun( @rdivide, total, t );

idl = skm( normTotal, 100, model, skmParams )==3;

% Remove rising/falling edges that may have low SNR and inaccurate FRET.
idl = imerode(idl, [1 1 1]);


idl(isnan(idl)) = 0;



end %function skmTotal
    
    