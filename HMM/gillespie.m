function [states,times] = gillespie(model, endTime)
% GILLESPIE  Stochastic simulation algorithm
%
%   [STATES,TIMES] = GILLESPIE(MODEL, ENDTIME) performs a stochastic simulation
%   of MODEL (QubModel object) using the Gillespie algorithm (direct method),
%   returning a sequence of STATES (1:model.nStates) and dwell times 
%   (TIMES, in milliseconds) that together are at least as long as ENDTIME.
%    
%   See also: simulate, simphotons, batchKinetics, QubModel.

%   Copyright 2018 Cornell University All Rights Reserved.



%% Pre-compute state time constants for Gillespie direct method
Q = model.rates;
nStates = size(Q,1);

Qtau    = zeros(1,nStates);  %mean dwell time for each state
Qcumsum = zeros(nStates);    %cumsum of probability of each possible exit from a state
for s=1:nStates
    Qtau(s)   = -1000 / sum( Q(s,:) );
    Qcumsum(s,:) = cumsum(  Q(s,:) ./ sum(Q(s,:))  );
end

% Sample initial state from initial probabilities distribution (p0)
p0 = model.p0/sum(model.p0);
curState = find( rand <= cumsum(p0), 1 );



%% Gillespie algorithm

% If available, use the compiled MEX version, which is 10-times faster.
try
    [states,times] = gillespiex(curState, endTime, Qtau, Qcumsum);
    return;
catch
end

%--- Otherwise use a fallback MATLAB implementation.
states = [];
times  = [];
cumTime = 0;

while cumTime<endTime,

    % Randomly sample the time until the next transition.
    dwellTime = Qtau(curState) .* log(rand);  %in ms

    % Randomly sample final state with probabilities calculated as the
    % fraction of all possible rate constants exiting current state.
    nextState = find( rand<=Qcumsum(curState,:), 1 );  %'first' is default

    states(end+1) = curState;
    times(end+1) = dwellTime;

    cumTime = cumTime + dwellTime;
    curState = nextState;

end %while not enough dwells to fill trace


