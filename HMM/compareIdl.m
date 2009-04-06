 function [accuracy,idl1,idl2] = compareIdl( dwtFilename1, dwtFilename2 )
% SKM  Crude model re-estimation using iterative idealization
% 
%   [DWT,NEW_MODEL] = SKM( DATA, SAMPLING, MODEL, params )
%   Optimizes the given FRET/kinetic MODEL soas to maximize the 
%   likelihood data given the model. DATA is a NxM matrix of
%   N FRET traces of M datapoints in length each. MODEL is a typical
%   model specification, as defined in createModel.m.  SKM
%   returns the optimal model (NEW_MODEL) and the idealization
%   with maximum likelihood (DWT).
%
%   NOTE that constraints on FRET values and stdev specified in
%   the model file WILL be enforced in the fitting with SKM.
%   Constraints on kinetics will be ignored!
%   
%   The following params may be specified: FIXME...
%    - maxItr (100):  maximum number of iterations before terminating
%    - convLL (1e-2): stop iterating when LL converges within this limit
%    - 
%   

% DEPENDS: idealize, forward_viterbi, countEvents??

% NOTE: aggregated states not yet supported!!!!

% NOTE for now we use the standard proceedure:
% 1. optimize for whole file and idealize
% 2. truncate each trace to last instance of non-zero FRET
% 3. optimize a model for each trace individually,
%    fixing mean FRET and stdev values.
% 4. idealize each trace with its optimal model...


% Get filenames of .dwt files to compare from user
if nargin==0,
    [f,p] = uigetfile('*.dwt','Select a DWT file');
    if f==0, return; end
    dwtFilename1 = [p f];
    
    [f,p] = uigetfile('*.dwt','Select another DWT file');
    if f==0, return; end
    dwtFilename2 = [p f];
else
    error('Too many input arguments');
end
    
% Load the .DWT files
[dwt1,sampling1,offsets1] = loadDWT(dwtFilename1);
[dwt2,sampling2,offsets2] = loadDWT(dwtFilename2);

if numel(offsets2)>1 && offsets2(2)==400,
    warning('applying a nasty hack.');
    offsets2 = offsets2*(25000/400);
end

% Convert DWTs to idealizations
idl1 = dwtToIdl( dwt1, offsets1 );
idl2 = dwtToIdl( dwt2, offsets2 );

% If files are of inconsistent sizes because they are taken
% at different framerates (for example when comparing simulated
% idealization to an estimate), expand the shorter one for to
% compensate.
if sampling1 ~= sampling2,
    % Determing the scaling factor
    sampFact = 1/(sampling2/sampling1);
    if sampFact<1,
        error('DWT2 must be bigger than DWT1 for scaling');
    end
    
    % Scale the smaller idealization
    idlFlat = idl1';
    idlFlat = idlFlat(:)';
    
    scaleIdx = repmat( 1:numel(idlFlat), [sampFact 1] );
    scaleIdx = scaleIdx(:);
    
    idlExpand = idlFlat( scaleIdx );
    idl1 = reshape( idlExpand, [size(idl2,2) size(idl2,1)] )';
end

% If the idealization in the first DWT ends before the second,
% truncate the second idl to the length of the first.
[nTraces,traceLen] = size(idl2);

for i=1:nTraces,
    e = find(idl1(i,:)>1,1,'last');
    if isempty(e),
        e=0;
    else
%         assert(~any(idl2(i,:)==0) );
        assert(~any(idl1(i,1:e)==0) );
    end
    idl2(i,e+1:end) = 0;
    idl1(i,e+1:end) = 0;
end

% Count differences
nDiff = sum( idl1(idl1~=0)~=idl2(idl2~=0) );
accuracy = 1-( nDiff/sum(idl1(:)~=0) );




 end







%%
% function idlFinal = dwtToIdl( dwt, traceLen )
% % NOTE how truncated data is handled?
% 
% nTraces = numel(dwt);
% idlFinal = zeros(traceLen,nTraces);
%     
% for dwtID=1:nTraces,
%     states = dwt{dwtID}(:,1);
%     times  = double(dwt{dwtID}(:,2));
%     assert(all(states>0));
% 
%     ends = [0; cumsum(times)];
%     for j=1:numel(states),
%         idlFinal( (ends(j)+1):ends(j+1),dwtID ) = states(j);
%     end
% end
% 
% idlFinal = idlFinal';
% 
% end


%%
function idlFinal = dwtToIdl( dwt, offsets )
% NOTE how truncated data is handled?

traceLen = median(diff(offsets));
idlFinal = zeros(offsets(end)+traceLen,1);
    
for dwtID=1:numel(dwt),
    states = dwt{dwtID}(:,1);
    times  = double(dwt{dwtID}(:,2));
    assert(all(states>0));

    ends = [0; cumsum(times)] +offsets(dwtID);
    for j=1:numel(states),
        idlFinal( ((ends(j)+1):ends(j+1)) ) = states(j);
    end
end

idlFinal = reshape(idlFinal,[traceLen,numel(idlFinal)/traceLen]);
idlFinal = idlFinal';

end














