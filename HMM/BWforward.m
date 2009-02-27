function [alphas,LLcarry] = BWforward( Bx, A, s0 ) 
%
% Usage: [alphas,LLcarry] = BWforward( O, A, mus, sigmas, s0 ) 
%
% Forward calculation of probabilities of observed sequence O up to time t and system is
% in particular state.  Probabilities divided into two terms, with relative probability
% given by alphas, and the bulk log-probability separated out at each time step in LLcarry.
% 
% O(t) is observed values, and A, B, s0 are given model paramters
%
% B is given by a gaussian with means mus and sigma sigmas
%
% DAB 2008.3.30

NT = size(Bx,1);
Nstates = size(A,1);
alphas = zeros(NT,Nstates);


% Setup initial conditions for forward calculations
% alphas(1,:) = s0 .* C .* exp( -((O(1)-mus).^2) ./ D );
alphas(1,:) = s0 .* Bx(1,:);
nrm = sum(alphas(1,:));
alphas(1,:) = alphas(1,:)/nrm;
LLcarry = log(nrm);

% Calculate emmission probabilities at each timepoint
% if the emmission came from each of the posssible states
% Bx = zeros( NT,Nstates );
% for i=1:Nstates
%     Bx(:,i) = C(i) .* exp(-  ((O-mus(i)).^2) ./ D(i)  );
% end

% Calculate alpha for each point in observation sequence
% alpha(t,i) = P( observations 1..t & state=i, given model[A,B,s0] )
for t = 2:NT
  % B(i) = P( observation t, given gaussian model[mu,sigma] in state=i )
%   B = C .* exp(-  ((O(t)-mus).^2) ./ D  );
  
  % alpha(time t, state j) = SUM_i( alpha(t-1,i)*A(i,j) ) * B(obs,j)
  % Probability of transitioning from each possible previous state i into
  % state j, using transition and observation probabilities in j and
  % collective probability of being in state i previously.
  % alpha*A is the same as al1*A(1,1) + al2*A(1,2) + al3*A(1,3)
  alphas(t,:) = (alphas(t-1,:)*A) .* Bx(t,:);
  nrm = sum(alphas(t,:));
  alphas(t,:) = alphas(t,:) / nrm;
  
  LLcarry = LLcarry + log(nrm);
end


