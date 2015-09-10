function [quality,randomScore] = weberQuality( base, registered, thresh )
% Calculates a Weber contrast score for random scampling of the target image.
% This provides a way to calculate the score expected for a random (poor)
% alignment. quality = weberScore/weberRandom(X,Y).
% Measure the "quality" of the alignment as the magnitude increase in
% score compared to a "random" alignment, which is approximated by
% scrambling the donor and acceptor field data. In practice, this gives
% similar scores to the mean of all scores from an alignment search across
% a wide parameter range. We do not just average the scores in the search
% because there may only be one parameter value or a small range. Quality
% scores above 1.1-1.15 (10-15% higher contrast than a random) are
% generally acceptable; any lower and the data may have other problems.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


N = 3; %number of repititions for averaging.
S = zeros(N,1);

for i=1:N,
    if i==1, %base case, no scrambling
        total = base(:)+registered(:);
    else
        total = base( randperm(numel(base)) ) + registered( randperm(numel(registered)) );
    end
    
    Ib = mean( total(total<thresh) );  %background intensity
    S(i) = ( mean(total(total>thresh)) -Ib ) / Ib;  %Weber contrast score
end

randomScore = mean(S(2:end));
quality = S(1) / randomScore;


end %FUNCTION weberRandom

