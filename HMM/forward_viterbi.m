function [vPath,vLL,tLL] = forward_viterbi(lsp, ltp, lep)
% FORWARD_VITERBI   Decodes sequence of states from observations (HMM)
%
%   [PATH,LL] = FORWARD_VITERBI(P0, A, B) finds the sequence of states PATH
%   that maximizes the probability of the data given the model (likelihood).
%   P0 is a Mx1 vector of the initial probabilities of each state.
%   A is the MxM transition probability matrix. A(1,2) is the probability at
%   each time step of transitioning from state 1 to state 2.
%   B is the MxN matrix of the probability at each time (cols) of observing the
%   data given the model, for each possible state (rows).
%   All inputs must be the natural logarithm of probabilities.
%
%   This functional automatically calls the compiled .mex version
%   forward_viterbix if availabile.
%
%   See also: idealize, skm.

%   Nice explainations:  http://www.comp.leeds.ac.uk/roger/HiddenMarkovModels/
%   html_dev/viterbi_algorithm

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Set a minimal value to probabilities to prevent underflow.
lsp(lsp==-Inf) = -1e10;
ltp(ltp==-Inf) = -1e10;
lep(lep==-Inf) = -1e10;


% Run the faster binary (mex) version if possible
try
    [vPath,vLL] = forward_viterbix( double(lsp),double(ltp),double(lep) );
    return;
    
catch err
    if strcmpi(err.identifier, 'MATLAB:UndefinedFunction'),
        disp('fowrard_viterbi: MEX not found. Falling back on the MATLAB version');
    else
        disp('WARNING: error in fowrard_viterbix.mex');
        disp(err);
    end
end


[nStates,nObs] = size( lep );

% Setup data structures used for finding the most likely path.
delta = zeros( nStates, nObs );  %partial probabilities
psi   = zeros( nStates, nObs );  %back pointers (most likely previous state)
% Maximal probability (and best path) that ends in state=i at time=t
% "If I am here, by what route is it most likely I arrived?"

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




