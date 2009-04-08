function [vPath,vLL,tLL] = forward_viterbi(start_p, trans_p, emit_p)
% FORWARD_VITERBI   Decodes sequence of states from observations (HMM)
%
%   [PATH,LL] = FORWARD_VITERBI(P0, A, B)
%   Finds the sequence of states (Viterbi PATH) with that is most
%   consistent with both the given observations and markov model
%   (high probability, LLv).  P0 are initial probabilities for each state,
%   A is the transition probability matrix, and B is a NxM matrix of
%   the probability the observation at time=t (1..N) was generated by
%   state s (1..M).

% Nice explainations:  http://www.comp.leeds.ac.uk/roger/HiddenMarkovModels/
% html_dev/viterbi_algorithm

% dbstop if error;

% Add a small number to all probabilities to prevent overflow (-Inf LL)
start_p = start_p+eps;
trans_p = trans_p+eps;
emit_p  = emit_p+eps;

%
[nStates,nObs] = size(emit_p);
assert( length(start_p)==nStates );

delta = zeros( nStates, nObs );  %partial probabilities
psi   = zeros( nStates, nObs );  %back pointers (most likely previous state)
% Maximal probability (and best path) that ends in state=i at time=t
% "If I am here, by what route is it most likely I arrived?"

lsp = log( start_p );
ltp = log( trans_p );
lep = log( emit_p  );

% Initiation
delta(:,1) = lsp + lep(:,1);

% Induction: calcualte probability at each timepoint
for t=2:nObs,
    %delta_t(j) = MAX_i[ delta_t-1(i) * A(i,j) ]   * B(j)
    
    for j=1:nStates  %current state
        
        % How likely am I to traverse all the way to step t-1,
        % then transition to j, and see the current observation?
        pCurr = delta(:,t-1) + ltp(:,j) + lep(j,t);
        % Which of the possible previous (source) states is most likely
        % given this information?
        [delta(j,t),psi(j,t)] = max(pCurr);
        
        %[valmax,argmax] = max(pCurr);
        % Most likely previous (source) state if current state is j
        %psi(j,t) = argmax;  
        % Prob. of being in state j at time t, given the optimal state
        % sequence up to this point.
        %delta(j,t) = valmax;  
    end
end

% Termination: find most likely endpoint
[valmax,argmax] = max(delta(:,end));
vLL = valmax/nObs;
tLL = sum(delta(:));

% Backtrace to determine most likely path to beginning
vPath = zeros(1,nObs);
vPath(nObs) = argmax;

for t=nObs-1:-1:1,
    vPath(t) = psi( vPath(t+1), t+1 );
end




