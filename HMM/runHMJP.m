function [idl,optModel,chain] = runHMJP(data, model, optsIn)
% runHMJP  Hidden Markov Jump Process model optimizer
%
%   [idl,optModel] = runHMJP(DATA, MODEL, OPTIONS)
%   Uses the Hidden Markov Jump Process approach to assign states and
%   optimize model parameters with DATA being a Traces object containing
%   the donor and acceptor fluoerscence traces as input and MODEL being a
%   QubModel object; only the number of states is used.
%   idl is the maximum a posteriori (MAP) output trajectory.
%   optModel is the output QubModel object with MAP estimate parameters.
%
%   The data should be curated so that there is no blinking (transitions
%   to non-fluorescent states) and the model should not include a dark
%   state. Frames after bleaching are removed with a threshold of 0.1.
%
%   OPTIONS is a struct with any additional settings (all fields optional):
%     .alpha    uniformization: determines further refinements of jump times within a frame period
%     .beta     higher number = slow rates. 1/(beta*eta)=peak escape rate in prior
%     .eta      gamma distribution shape parameter: 4=peaked prior, 2=exp prior.
%     .HMC_eps  Hamiltonian Monte Carlo integration step size
%     .HMC_L    Hamiltonian Monte Carlo number of Leap-frog integration steps.n.
%     .removeBleaching: remove final dwell in dark state (true).
%
%   Algorithm details can be found in the following publications:
%     "Generalizing HMMs to Continuous Time for Fast Kinetics: Hidden Markov Jump Processes"
%       https://doi.org/10.1101/2020.07.28.225052
%     "Extraction of rapid kinetics from smFRET measurements using integrative detectors"
%       https://doi.org/10.1101/2020.08.28.267468
%
%   See also: batchKinetics, runParamOptimizer, chainer_main.

%   Copyright 2022 St Jude Children's Research Hospital.

% This is a very preliminary implementation. Future steps:
% 1. return continuous time trajectories (MAP per trace).
% 2. repeat until convergence instead of fixed number of iterations.
% 3. return/display distributions of parameters across traces.
% 4. return/display probability distributions??
% 5. Parallelize; eventually allow submission to a cluster.
% 6. Incorporate stroboscopic info (if any).


%% Process input arguments
narginchk(2,3);

if ischar(data), data=loadTraces(data); end
assert( isa(data,'Traces'), 'First argument must be Traces object or path to .traces file' );

if ischar(model), model=QubModel(model); end
assert( isa(model,'QubModel'), 'First argument must be QubModel object or path to .model file' );


%% Configure HMJP sampler

opts = hmmopt(mfilename);
opts.units.space = 'nm';
opts.units.Int   = 'photons';
opts.units.time  = 's';
opts.M           = model.nStates;
opts.mu_back_D   = 0;  % Background fluorescence levels
opts.mu_back_A   = 0;

opts = mergestruct( opts, optsIn );


%% Run HMJP, requesting maximum likelihood point estimates
idl = zeros( data.nTraces, data.nFrames );
optModel(data.nTraces) = QubModel;
dt = data.sampling/1000;  %frame interval in seconds.

for traceID=1:data.nTraces
    
    if opts.exclude(traceID), continue; end
    
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
    
    %while not converged...

        % Run HMJP sampler for chain_length iterations
        chain = chainer_main(chain, opts.maxItr, [],true,true);

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
    
    % Sort states by increasing mean FRET value to fix label switching
    mu = chain.mu_A(:,best) ./ (chain.mu_A(:,best)+chain.mu_D(:,best));
    [mu,sortIdx] = sort(mu);
    optModel(traceID) = copy(model);
    optModel(traceID).p0 = chain.IT(best,sortIdx);
    optModel(traceID).mu = mu;
    %optModel(traceID).sigma = ???
    optModel(traceID).rates = chain.Q(sortIdx,sortIdx,best);
    
    % Discretize trajectory and construct dwell-time matrix.
    starts = floor(chain.t_t{best}/dt)+1;  %convert frame numbers
    [starts,idx] = unique(starts);
    states = sortIdx(chain.t_s{best}(idx));
    
    for j=1:numel(starts)-1
        idl(traceID, starts(j):starts(j+1)) = states(j);
    end
    idl(traceID, starts(end):floor(opts.T_f/dt)) = states(end);
    
end  %for each trace



end %function runHMJP



