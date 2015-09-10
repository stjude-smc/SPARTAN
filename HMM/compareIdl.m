 function [accuracy,idl1,idl2] = compareIdl( dwtFilename1, dwtFilename2 )
% compareIdl    Calculates percent similarity of idealizations.
% 
%   [ACC] = compareIdl( DWT1, DWT2 )
%   Calculates the percentage (ACC) of FRET datapoints that is the same
%   in the idealizations saved in DWT1 and DWT2. DWT1 and DWT2 are filnames
%   of QuB-format dwell-time files (.dwt). The files may have different
%   framerates, as long as DWT2 has a faster framerate than DWT1. This
%   feature is useful for comparing the results of idealization algorithms
%   with the simulated state sequence (true result).

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% DEPENDS: loadDWT, dwtToIdl.
% TODO: .


% Get filenames of .dwt files to compare from user
if nargin<1,
    [f,p] = uigetfile('*.dwt','Select a DWT file');
    if p==0, return; end
    dwtFilename1 = [p f];
end

if nargin<2
    [f,p] = uigetfile('*.dwt','Select another DWT file');
    if p==0, return; end
    dwtFilename2 = [p f];
end

assert(  exists(dwtFilename1,'file') && exists(dwtFilename2,'file')  );

% Load the .DWT files
[dwt1,sampling1,offsets1] = loadDWT(dwtFilename1);
[dwt2,sampling2,offsets2] = loadDWT(dwtFilename2);

% Convert DWTs to idealizations
traceLen = median(diff(offsets1));
idl1 = dwtToIdl( dwt1, traceLen, offsets1 );
traceLen = median(diff(offsets2));
idl2 = dwtToIdl( dwt2, traceLen, offsets2 );

% If files are of inconsistent sizes because they are taken
% at different framerates (for example when comparing simulated
% idealization to an estimate), expand the shorter one to compensate.
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












