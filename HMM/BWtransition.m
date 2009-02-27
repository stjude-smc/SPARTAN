function [eps,alphas,LLcarry] = BWtransition( O, A, mus, sigmas, p0 )
%
% Usage: [eps,alphas,LLcarry] = BWtransition( O, A, mus, sigmas, p0 )
%
% Calculates probability of transition between states for each time
% point given model [A,mus,sigmas,p0] following adapted Baum-Welch algorithm
%
% DAB 2008.3.30

NT = length(O);
Nstates = size(A,1);

% Calculate emmission probabilities at each timepoint
% if the emmission came from each of the posssible states
% C = 0.3989422804014327./sigmas; %gaussian leading coefficient
% C = 1./sqrt(2*pi*sigmas); %gaussian leading coefficient
% D = (2*sigmas.^2);  %exponential denominator

Bx = zeros( NT,Nstates );
for i=1:Nstates
%     Bx(:,i) = C(i) .* exp(-  ((O-mus(i)).^2) ./ D(i)  );
    Bx(:,i) = normpdf( O, mus(i), sigmas(i) )/6;
end

% Calculate forward (alpha) and backward (beta) probabilities
% forward  = prob. of obs. up to time t and the current state is i
% backward = prob. of obs. following time t and the current state is i
[alphas LLcarry] = BWforward(Bx,A,p0);
[betas LL2] = BWbackward(Bx,A);


% Don't try to optimize this loop
% Epsilon = prob. of getting to state i (alpha) with state j following (beta)
%   at time t, with the transition between these states and the emission
%   at that state (prob. of observation) included.
%   probability of being in state i at time t and going to state j,
%   given the current model 
% 
% In other words, Epsilon is a the transition probability matrix at time=t
%   given complete knowledge of the sequence preceeding and following.

eps = cell(NT-1,1);
es = zeros(Nstates,Nstates);
for t = 1:NT-1  %for each datapoint
  for i = 1:Nstates
    for j = 1:Nstates  %each possible start/end state pair
      es(i,j) = alphas(t,i)*betas(t+1,j) * A(i,j) .* Bx(t+1,j);
    end
  end
  eps{t} = es / sum(es(:)); %normalize
end

% eps should be a 3D matrix instead of cell array for speed!
