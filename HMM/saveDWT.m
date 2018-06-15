function saveDWT(filename, varargin)
% SAVEDWT  Saves a QuB-format idealization dwell time file
%     
%   saveDWT( FILE, DWT, OFFSETS, MODEL, SAMPLING )
%   Saves the Dwell-Times data (DWT) to FILENAME.  SAMPLING is in ms.
%   MODEL is a Nx2 matrix of mean and stdev values for each of N  states.
%   DWT is a cell array containing dwell-time matrices on N traces.
%   Each cell element is a 2xM matrix, with state+dwell time pairs in
%   each row (see saveDWT.m).  Times are in *frames* and states are 1-based.
%   All necessary conversion are made within this function.
%   OFFSETS is an array of indexes for the start of each segment
%   into the parent data file (.qub.txt).
%
%   saveDWT( FILE, IDL, MODEL, SAMPLING )
%   Save a .dwt file using the state assignment matrix IDL instead of a cell
%   array of dwell-time series.
%
%   See also: loadDWT.

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% Example of the file format used:
% 
% Segment: {N} Dwells: {M} Sampling(ms): {R} Start(ms): {I} ClassCount: {J} Mu1 Std1 Mu2 Std2
% Segment: 1 Dwells: 15 Sampling(ms): 40 Start(ms): 0 ClassCount: 4 0.010000 0.050000 0.240000 0.061000 0.390000 0.061000 0.560000 0.061000
% 3	2800
% 1	320
% 0	40
% 1	400
% 3	600
% 1	40


narginchk(4,5);
nargoutchk(0,0);


% Process dwell-time series input arguments
if iscell(varargin{1})
    [dwt, offsets, model, sampling] = varargin{:};

elseif isnumeric(varargin{1})
    [idl, model, sampling] = varargin{:};
    [dwt,offsets] = idlToDwt(idl);
else
    error('Invalid idealization input. Should be cell array or assignment matrix')
end


% If model is a qub model structure, convert it into the expected form.
if isa(model,'QubModel') || isstruct(model)
    model = [to_col(model.mu) to_col(model.sigma)];
end


% verify input arguments
assert( numel(dwt)==numel(offsets), 'Offsets do not match idealization' );


% Switch to row order so that (:) creates a sequence of mean+stdev pairs
if ~iscell(model),
    assert( size(model,2)==2, 'Invalid model array shape' );
    
    fretValues = model(:,1);
    if any( diff(fretValues)<=0 ),
        warning('saveDWT: FRET values are non-increasing!');
    end
end
m = model;


% Save idealization information to file
fid = fopen(filename,'w');
disp( ['Saving to ' filename] );

for segId=1:numel(dwt),
    
    % Get the model used for this trace, if applicable.
    if iscell(model),
        m = model{segId};
    end
    nStates = size(m,1);
    
    % Get dwell-time for the current segment, convert to ms.
    segment = dwt{segId};
    offset  = offsets(segId)*sampling;
    segment(:,2) = segment(:,2)*sampling;
    
    % Convert to zero-based class numbers used by the file format.
    segment(:,1) = segment(:,1)-1;
    assert( all(segment(:,1)<nStates), 'Invalid state number' );
    
    % Write segment data to file. (sprintf is faster than fprintf)
    text = [ sprintf('Segment: %d Dwells: %d Sampling(ms): %d Start(ms): %d ClassCount: %d', ...
                 segId, size(segment,1), sampling, offset, nStates) ...;
             sprintf(' %f', m' ) '\n' sprintf('%d\t%d\n', segment') ];
    fprintf(fid,text);
end

fclose(fid);

end  % function LoadDWT
