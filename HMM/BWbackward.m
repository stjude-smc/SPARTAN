function [betas,LLcarry] = BWbackward( Bx, A ) 
%
% Usage: [betas,LLcarry] = BWbackward( O, A, mus, sigmas ) 
% 
% Backward calculation of probality of observed sequence and state (see BWforward for more info)
%
% O(t) is observed values, and A, B, s0 are given model paramters
%
% B is given by a gaussian with means mus and sigma sigmas
%
% DAB 2008.3.30

NT = size(Bx,1);
Nstates = size(A,1);

betas = zeros(NT,Nstates);
betas(NT,:) = 1;

LLcarry = 0;


% Calculate backward probabilities
% beta(t,i) = P( observations t+1..end, given state=i & model[A,B,s0] )
for t = NT-1:-1:1
  % B(i) = P( observation t, given gaussian model[mu,sigma] in state=i )
%   B = C .* exp(-  ((O(t+1)-mus).^2) ./ D  );
  
  betas(t,:) = A* (betas(t+1,:).*Bx(t+1,:))';
%   betas(t,:) = (betas(t+1,:).*B) * A;
  nrm = sum(betas(t,:));
  betas(t,:) = betas(t,:)/nrm;
  
  LLcarry = LLcarry + log(nrm);
end
