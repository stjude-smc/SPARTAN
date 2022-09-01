function [idl,optModel,chain] = runHMJP(data, model, optionsInput)
% runHMJP  Hidden Markov Jump Process model optimizer
%
%   [idl,optModel] = runHMJP(DATA, MODEL, OPTIONS)
%   Uses the Hidden Markov Jump Process approach to assign states and
%   optimize model parameters with DATA being a Traces object containing
%   the donor and acceptor fluoerscence traces as input and MODEL being a
%   QubModel object with the starting parameter values.
%   idl is the maximum a posteriori (MAP) output trajectory.
%   optModel is the output QubModel object with MAP estimate parameters.
%
%   OPTIONS is a struct with any additional settings (all fields optional):
%     .maxIter:  maximum number of iterations (150)
%     .convLL:   termination tolerance in LL values.
%     .convGrad: termination tolerance for parameter step size.
%     .verbose:  print information about each iteration for debugging.
%     .updateModel: modify input model object in each iteration (false).
%                Enables rate constants to be viewed during optimization.
%     .removeBleaching: remove final dwell in dark state (true).
%
%   Algorithm details can be found int he following publications:
%  
%   "Generalizing HMMs to Continuous Time for Fast Kinetics: Hidden Markov Jump Processes"
%   https://doi.org/10.1101/2020.07.28.225052
%
%   "Extraction of rapid kinetics from smFRET measurements using integrative detectors"
%   https://doi.org/10.1101/2020.08.28.267468
%
%   See also: batchKinetics, runParamOptimizer, chainer_main.

%   Copyright 2022 St Jude Children's Research Hospital.

narginchk(2,3);
[idl,optModel,chain] = deal([]);

% verbose = options.verbose;  %show intermediate results from HMJP.


%% Configure HMJP sampler

opts.units.space = 'nm';
opts.units.Int   = 'photons';
opts.units.time  = 's';

chain_length = 500;

% Hyperparameters
opts.M        = 2;  %model.nStates;
opts.alpha    = 1.01;      % Uniformization; determines further refinements of jump times within a frame period
opts.beta     = 5;      %higher number = slow rates. 1/(beta*eta)=peak escape rate in prior
opts.eta      = 4;      %4=peaked prior, 2=exp prior.
opts.HMC_eps  = 0.01;   % Hamiltonian Monte Carlo integration step size
opts.HMC_L    = 50;     % Hamiltonian Monte Carlo number of Leap-frog integration steps.

% Background fluorescence levels
opts.mu_back_D   = 0;
opts.mu_back_A   = 0;



%% Run HMJP, requesting maximum likelihood point estimates
% dwt = cell( data.nTraces, 1 );
idl = zeros( data.nTraces, data.nFrames );
dt = data.sampling/1000;  %frame interval in seconds.

for traceID=1:data.nTraces
    
    % Trace data preprocessing:
    % remove photobleached portion of trace and
    % clip negative values that cause the algorithm to crash.
    % FIXME: for now ignore any stroboscopic information.
    % FIXME: need to segment to avoid blinking events.
    idxEnd = find( data.fret(traceID,:)>0.1, 1, 'last' )-1;
    opts.Int_D    = max(0, data.donor(traceID,1:idxEnd) );
    opts.Int_A    = max(0, data.acceptor(traceID,1:idxEnd) );
    
    opts.t_left   = data.time(1:idxEnd)/1000;  %left edge of each frame time
    opts.t_right  = data.time(2:idxEnd)/1000;    %right edge of each frame time
    opts.t_right  = [opts.t_right opts.t_right(end)+dt];
    opts.T_i      = opts.t_left(1);         %start point in seconds
    opts.T_f      = opts.t_right(end);      %end point in seconds
    
    opts.phi_D    = 1; % mean of the mean_st DONOR
    opts.psi_D    = sum(opts.Int_D)/(2*(opts.T_f-opts.T_i));
    opts.phi_A    = 1; % mean of the mean_st ACCEPTOR
    opts.psi_A    = sum(opts.Int_A)/(2*(opts.T_f-opts.T_i));
    opts.alpha_prop_D  = opts.psi_D;
    opts.alpha_prop_A  = opts.psi_A;
    opts.theta_prop_D  = opts.psi_D;
    opts.omega_prop_A  = opts.psi_A;
    
    % Initialize HMJP sampler
    chain = chainer_main( [], 0, opts, true, [] );
    
    %for j=1:100 %maxIterations

        % Run HMJP sampler for chain_length iterations
        chain = chainer_main(chain, chain_length, [],true,true);

        % Truncate chain if too large?
        %if chain.sizeGB > 1.0
        %    chain = chainer_main(chain,-fix(chain.length/2),[],true,[]);
        %end

        % Plot current results from current iteration.
        % These functions don't seem to be included in Zeliha's folder.
        %export_chain(chain, 1 , ['tag_mean_sim_NEW_',num2str(tag),'_',num2str(tag_case),'_',num2str(length(chain.i))]);
        %plot_escrate;
        %plot_diag_prob;
        %plot_trajectories;
        %plot_state_level;
        %drawnow;
        
        % Test for convergence here
        %if abs(mean(chainer.MAP(1:end))) - abs(mean(chainer.MAP(1:end-100))) < tol
        %    run additional 3000-4000 steps...
        %    break;
    %end
    
    % Choose maximum a posteriori iteration and save result
    [~,best] = min(chain.MAP);
    
    % Fix label switching: sort by increasing model mean FRET value
    mu = chain.mu_A(:,best) ./ (chain.mu_A(:,best)+chain.mu_D(:,best));
    [mu,sortIdx] = sort(mu);
    
    % FIXME: hacky way of dealing with the dark state...
%     optModel.rates(2:end,2:end) = chain.Q(sortIdx,sortIdx,best);
%     optModel.p0    = [0 chain.IT(best,sortIdx)];
%     optModel.mu    = [0; mu];
    
    % Discretize trajectory and construct dwell-time matrix.
    % FIXME: assumes trace begins at time=0!
    starts = floor(chain.t_t{best}/dt)+1;  %convert frame numbers
    [starts,idx] = unique(starts);
    states = sortIdx(chain.t_s{best}(idx))  +1;  %+1 for dark state!
    
    for j=1:numel(starts)-1
        idl(traceID, starts(j):starts(j+1)) = states(j);
    end
    idl(traceID, starts(end):floor(opts.T_f/dt)) = states(j);
    
    % Deal with inverted states
    
    
end  %for each trace

% For now, don't try to report rates etc.
% This should report one model per trace like MIL Separately!
optModel = copy(model);


end %function runHMJP

