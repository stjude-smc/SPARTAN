function [Etot,gamma,LL,alpha,beta] = BWtransition( data, A, mu, sigma, p0 )
% Forward-backward algorithm for hidden Markov modeling
%
%   [eps,alphas,LL] = BWtransition( DATA, A, mu, sigma, p0 )
%   DATA is a vector of experimental data values (FRET trace).
%   A is the transition probability matrix.
%   MU/SIGMA are the mean/stdev of each state's Gaussian emission distribution.
%   P0 are the initial probabilities of each state.
%
%   See also: batchKinetics, bwOptimize, forwardBackward.

%   Copyright 2008-2018 Cornell University All Rights Reserved.


narginchk(5,5);
nargoutchk(1,5);
p0 = reshape(p0, 1, numel(p0));  %ensure row vector


%% Calculate emmission probabilities at each timepoint
% The division by 6 keeps the observation probabilities <1 (negative LL).
nFrames = numel(data);
nStates = size(A,1);

Bx = zeros(nFrames, nStates);
for i=1:nStates
    Bx(:,i) = normpdf( data, mu(i), sigma(i) );
end


% Run compiled version if available. Fallback on Matlab otherwise.
try
    [LL,alpha,beta,gamma,Etot] = forwardBackward(p0, A, Bx);
    return;
catch
end


%% Forward probabilities
% alpha(t,i) = P( observations 1..t & state(t)=i | model )
alpha = zeros(nFrames, nStates);
nrm = zeros(nFrames,1);

alpha(1,:) = p0 .* Bx(1,:);
nrm(1) = 1./sum(alpha(1,:));
alpha(1,:) = alpha(1,:) * nrm(1);

for t=2:nFrames
    %alpha(t,j) = SUM_i( alpha(t-1,i)*A(i,j) ) * B(t,j)
    alpha(t,:) = (alpha(t-1,:) * A) .* Bx(t,:);

    nrm(t) = 1./sum(alpha(t,:));  %normalization prevents float underflow
    alpha(t,:) = alpha(t,:) * nrm(t);
end

LL = -sum( log(nrm) );


%% Backward probabilities
% beta(t,i) = P( observations t+1..end | state(t)=i & model )
beta = zeros(nFrames, nStates);
beta(nFrames,:) = 1;

for t = nFrames-1:-1:1
    %beta(t,i) = SUM_j(  A(i,j) * B(t+1,j) * beta(t+1,j)  )
    beta(t,:) = A *  ( beta(t+1,:) .* Bx(t+1,:) )'  * nrm(t);
end


%% Calculate probability of being in each state at each time:
% gamma(t,i) = P( state(t)=i | all obseravtions & model )
% SUM_t(gamma) is the expected number of times each state is occupied.
% Summing E terms below is equivalent, but slower.
gamma = alpha .* beta;
gamma = bsxfun( @rdivide, gamma, sum(gamma,2) );  %normalized at each time t


%% E = Joint prob. of being in state i at time t AND state j and time t+1.
% In other words, the transition probability t->t+1 given all observed data.
% SUM_t(E) is the expected number of each type of transition.
% See eq. 37 of Rabiner 1989.
Etot = zeros(nStates);
es = zeros(nStates);

for t=1:nFrames-1  %for each datapoint
    %E(t,i,j) = alpha(t,i) * A(i,j) * B(t+1,j) * beta(t+1,j) / norm(E(i,j))
    for i=1:nStates
        es(i,:) = alpha(t,i) * A(i,:) .* Bx(t+1,:) .* beta(t+1,:);
    end
    Etot = Etot + es / sum(es(:)); %normalize
end


