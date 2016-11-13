function [bestIdl,bestModel] = runVbFret(data) %, vb_opts_in, PriorPar_in)
% runVbFret   Idealize traces using the vbFRET program.
% 
%   [IDL,MODEL] = runVbFret(DATA) run's Bronson et al's vbFRET program to
%   find optimal models and state assignment traces for each input FRET trace.
%   DATA may be a matrix, TracesFret object, or .traces filename.
%   IDL is a state assignment matrix, same size as DATA.fret.
%   MODEL is a cell array of structs with optimal vbFRET model parameters.
%   
%   Modified from the vbFRET_no_gui.m function provided with vbFRET.
%   For algorithm information, see:  http://vbfret.sourceforge.net/
%
%   See also: batchKinetics, bwOptimize, skm, idealize.

%   Copyright 2016 Cornell University All Rights Reserved.

narginchk(1,1);

if ischar(data),
    data = loadTraces(data);
end
if isa(data,'TracesFret')
    fret = data.fret;
elseif isnumeric(data)
    fret = data;
else
    error('Invalid input data.');
end



%%%%%%%%%%%%%%%%%%%%%
% parameter settings
%%%%%%%%%%%%%%%%%%%%%

D = 1;  % analyze data in 1D (FRET, not fluorescence)
minStates = 2;  %minimum number of states to try fitting
maxStates = 5;  %maximum number of states to try fitting
maxRestarts = 10;  %maximum number of restarts 

% Priors for optimization algorithm
PriorPar.upi  = 1;
PriorPar.mu   = .5*ones(D,1);
PriorPar.beta = 0.25;
PriorPar.W    = 50*eye(D);
PriorPar.v    = 5.0;
PriorPar.ua   = 1.0;
PriorPar.uad  = 0;
%PriorPar.Wa_init = true;

% Optimization settings
vb_opts.maxIter     = 100;   %stop after vb_opts iterations if program has not yet converged
vb_opts.threshold   = 1e-5;  %should this be a function of the size of the data set??
vb_opts.displayFig  = 0;     %display graphical analysis
vb_opts.displayNrg  = 0;     %display nrg after each iteration of forward-back
vb_opts.displayIter = 0;     %display iteration number after each round of forward-back
vb_opts.DisplayItersToConverge = 0;  % display number of steps needed for convergance

% FIXME: use mergestruct for defaults above


nTraces = size(fret,1);

% bestMix = cell(nTraces,maxStates);
bestOut  = cell(nTraces,maxStates);
outF     = -inf*ones(nTraces,maxStates);  %best score for each condition
% best_idx = zeros(nTraces,maxStates);



%%%%%%%%%%%%%%%%%%%%%%%%
%run the VBEM algorithm
%%%%%%%%%%%%%%%%%%%%%%%%
tic;

% Iterate over traces
for n=1:1  %nTraces,
    % Truncate trace for faster processing
    trace = fret(n,:);
    trace = trace(trace~=0);  %must be a row vector
    
    % For each total number of states to try
    parfor k=minStates:maxStates
        init_mu = (1:k)'/(k+1);  %uniform distribution for first try
        i=1;  %restart number
        maxLP = -Inf;  %max score across restarts in this condition
        
        % Repeat maxRestarts times, fitting with random state means.
        while i<maxRestarts+1
            if k==1 && i>3, break; end  %??

            if i>1, init_mu=rand(k,1); end   %Random starting emission means.

%             fprintf('Working on inference for restart %d, k%d of trace %d...\n',i,k,n);

            % Initialize gaussians
            % Note: try-catch needed for the K-means algorithm. If it doesn't
            % initialze, vbFRET_VBEM may crash.
            try
                % Construct Gaussian mixture emission models
                [mix] = get_mix(trace',init_mu);
                
                % Run optimization algorithm (vbFRET)
                [out] = vbFRET_VBEM(trace, mix, PriorPar, vb_opts);
                
            catch e
                disp('There was an error, repeating restart.');
                disp(e.message);
                continue;
            end
            
            % Only save the iterations with the best out.
            % FIXME: this can be simplified by only keeping the optimal model.
            if out.F(end) > maxLP
                maxLP = out.F(end);
                %bestMix{n,k} = mix;
                bestOut{n,k} = out;
                outF(n,k)=out.F(end);
                %best_idx(n,k) = i;
            end
            
            i=i+1;
            
        end %for each restart
        
    end %for each number of states
    
end %for each trace



%%%%%%%%%%%%%%%%%%%%%%%% VBHMM postprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze accuracy and save analysis
disp('Idealizing using optimal models...')


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get idealized data fits
%%%%%%%%%%%%%%%%%%%%%%%%%%
% z_hat = cell(nTraces,maxStates);  %state assignments
% x_hat = cell(nTraces,maxStates);  %idealized fret trace
% 
% for n=1:nTraces
%     for k=minStates:maxStates
%         %[z_hat{n,k},x_hat{n,k}] = chmmViterbi(bestOut{n,k},FRET{n}(:)');
%         z_hat{n,k} = chmmViterbi( bestOut{n,k}, data.fret(n,:) );
%     end
% end



% Use the highest-scoring model to idealize each trace.
% FIXME: model may not have a zero state and may misassign points after
% bleaching.
bestModel = cell(nTraces,1);
% bestIdl   = cell(nTraces,1);
bestIdl = zeros(size(fret));

for n=1:1  %nTraces,
    idxBest = find( outF(n,:)==max(outF(n,:)), 1,'first' );
    bestModel{n} = bestOut{n,idxBest};
    
    bestIdl(n,:) = chmmViterbi( bestModel{n}, fret(n,:) );  %second parameter is idealized trace
end

disp(toc);

% Description of outputs:
% States are not necessarily in increasing order of FRET values.
% 
% Wa   = unnormalized transition probability matrix (A).
% Wpi  = unnormalized initial probabilities (pi).
% beta = internal (not useful for us).
% m    = mean FRET state values.
% W,v  = parameters used to construct the covariance matrix.
% F    = model likelihood score (array of each iteration).
%
% FIXME: need to output in an expected form...


disp('...done.')    

end %function runVbFret

