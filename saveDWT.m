function saveDWT(filename, idealization, offsets, model, sampling )
% SAVEDWT  Saves a QuB-format idealization dwell time file
%     
%   SaveDWT( FILE, DWT, OFFSETS, MODEL, SAMPLING )
%   Saves the Dwell-Times data (DWT) to FILENAME.  SAMPLING is in ms.
%   MODEL is a Nx2 matrix of mean and stdev values for each of N  states.
%   DWT is a cell array containing dwell-time matrices on N traces.
%   Each cell element is a 2xM matrix, with state+dwell time pairs in
%   each row (see saveDWT.m).  Times are in *frames* and states are 1-based.
%   All necessary conversion are made within this function.
%   OFFSETS is an array of indexes for the start of each segment
%   into the parent data file (.qub.txt).
%
%   All arguments are required!

% Example datafile
% Segment: {N} Dwells: {M} Sampling(ms): {R} Start(ms): {I} ClassCount: {J}
% Mu1 Std1 Mu2 Std2 ...
% Segment: 1 Dwells: 15 Sampling(ms): 40 Start(ms): 0 ClassCount: 4 0.010000 0.050000 0.240000 0.061000 0.390000 0.061000 0.560000 0.061000
% 3	2800
% 1	320
% 0	40
% 1	400
% 3	600
% 1	40


% If model is a qub model structure, convert it into the expected form.
if isstruct(model),
    mu = model.mu(model.class);
    sigma = model.sigma(model.class);
    
    model = [ reshape(mu,numel(mu),1) reshape(sigma,numel(sigma),1) ];
end


% verify input arguments
assert( iscell(idealization), 'Idealization must be a cell array' );
assert( numel(idealization)==numel(offsets), 'Offsets do not match idealization' );
assert( size(model,2)==2, 'Incorrect model dimensions' );

nSegments = numel(idealization);


% Switch to row order so that (:) creates a sequence of mean+stdev pairs
if ~iscell(model),
    assert( size(model,2)==2, 'Invalid model array shape' );
    model = model'; % FRET in first row, stdev in second row.
    nStates = numel(model)/2;
    
    fretValues = model(1,:);
    if any( diff(fretValues)<=0 ),
        warning('saveDWT: FRET values are non-increasing!');
    end
else
    nStates = numel(model{1})/2;
end


% Save idealization information to file
fid = fopen(filename,'w');
disp( ['Saving to ' filename] );

for ctrSeg=1:nSegments,
    
    segment = idealization{ctrSeg};
    nDwells = size(segment,1);
    offset  = offsets(ctrSeg);
    
    % Substract 1 from state numbers (QuB uses 0-based)
    segment(:,1) = segment(:,1)-1;
    assert( all(segment(:,1)<nStates), 'Invalid state number' );
    
    % Convert frames to ms
    segment(:,2) = segment(:,2)*sampling;
    
    % Write segment header
    fprintf(fid, 'Segment: %d Dwells: %d Sampling(ms): %d Start(ms): %d ClassCount: %d', ...
                 ctrSeg, nDwells, sampling, offset*sampling, nStates );
    
    if iscell(model),
        m = model{ctrSeg}';
        m = m(:);
    else
        m = model(:);
    end    
    fprintf(fid, ' %f', m );
    fprintf(fid, '\n');
    
    % Write sequence of dwells
    fprintf(fid, '%d\t%d\n', segment');
end

fclose(fid);

end  % function LoadDWT
