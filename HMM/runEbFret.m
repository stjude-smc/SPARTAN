function [idlTotal,optModels,LL,selfanalysis] = runEbFret(data, dt, model, params)
% runEbFret  Empirical Bayes model optimization and idealization.
%   
%   [IDL,MODEL] = runVbFret(DATA,STATES) finds the most likely model given the
%   data using the empirical Bayes method (see below) and creates idealized
%   FRET traces (state assignment) using the Viterbi algorithm.
%   DATA is a TracesFret object or .traces filename.
%   STATES is a list of number of states (model complexity) to try.
%
%   "ebFRET is the successor of vbFRET, which it improves upon by performing
%   combined analysis of multiple time series at once. This approach, known as
%   empirical Bayes, is both more accurate and statistically robust. It also 
%   enables more advanced analysis use cases such as sub-population detection."
%   
%   This function requires the ebFRET package to be on the MATLAB path.
%   ebFRET can be downloaded here: https://ebfret.github.io/
%   This file is modified from ebFRET's run_ebayes.m and run_vbayes.m files.
%   
%   See also: batchKinetics, skm, idealize, bwOptimize.

%   Copyright 2016 Cornell University All Rights Reserved.

narginchk(2,Inf);

    

    wbh = waitbar(0,'Running ebFRET...');
    dL = 0.3;

%     ip = inputParser();
%     ip.StructExpand = true;
%     ip.addParamValue('analysis', [], @isnumeric);       
%     ip.addParamValue('series', [], @isnumeric);       
%     ip.addParamValue('threshold', 1e-4, @isnumeric);       
%     ip.addParamValue('max_iter', 100, @isscalar);       
%     ip.parse(varargin{:});
%     args = ip.Results;

    % Process input arguments
    if ischar(data), data=loadTraces(data); end
    if isa(data,'TracesFret'), data=data.fret; end
    if ~isnumeric(data), error('Invalid FRET data input'); end
    
    origSize = size(data);
    data = data(~params.exclude,:);
    
    model = copy(model); %we don't want to modify in place.
    nStates = model.nStates; %2:3;
    nTraces = size(data,1);

    % Remove outlier values
    data = min( max(data,-0.2), 1.2);
    
    % Truncate traces when donor blinks or bleaches.
    % Stretches of zeros (donor blinking) tend to cause vbayes to crash.
    % FIXME: run viterbi outside run_vbayes to fix this.
    xall = cell(nTraces,1);
    exclude = false(nTraces,1);
    
    for i=1:nTraces,
        lt = find(data(i,:)~=0,1,'last');
        xall{i} = double( data(i,1:lt)' );
        
        % Add zero-FRET noise so all traces have a zero-FRET state.
        xall{i} = [ xall{i}; randn(20,1)*0.05 ];
    end
    nFrames = cellfun(@numel,xall);
    
    % settings
    threshold = 1e-5;
    num_restarts = 2;
    max_iter = 100;

    % create analysis struct (priors) with default values.
    selfanalysis = init_analysis2(xall, nStates);

    % loop over different analysis sets
    for a = nStates
        
        % start empirical bayes iterations
        it = 1;
        while true
            % do variational bayes updates for posterior
            selfanalysis = run_vbayes(xall, selfanalysis, a, ...
                'restarts',  num_restarts, ...
                'threshold', threshold, ...
                'max_iter',  max_iter);
            
            % Exclude any problematic traces from further analysis.
            % These prevent convergence of ensemble parameters.
            exclude = selfanalysis(a).lowerbound./nFrames>3;  %FIXME: need a better threshold
            lb = selfanalysis(a).lowerbound(~exclude) ./ nFrames(~exclude) / sum(~exclude);
            L(it) = sum(lb);
            
            if it == 1
                %fprintf('it %02d   L %.5e\n', it, L(it));
            else
                dL = (L(it)-L(it-1)) / abs(L(it));  %relative change in likelihood score
                %fprintf('it %02d   L %.5e    dL %.2e\n', it, L(it), dL);
            end
            if ~ishandle(wbh) || ~isvalid(wbh)
                error('spartan:op_cancelled','Operation cancelled by user');
            end
            progress = max( it/max_iter, log10(dL)/log10(threshold) );
            waitbar(progress, wbh);

            % check termination crtieria
            if it>max_iter, break; end
            if it>1 && dL<threshold, break; end
            
            % run iterative empirical bayes update
            % (assuming constant posterior statistics)
            u = selfanalysis(a).prior;
            w = selfanalysis(a).posterior(~exclude);
            E = selfanalysis(a).expect(~exclude);
            selfanalysis(a).prior = ebfret.analysis.hmm.h_step(w, u, 'expect', E);

            % increment iteration counter
            it = it + 1;
        end
    end
    
    if any(exclude),
        fprintf('Excluded %d traces (%.0f%%) from analysis:\n', sum(exclude), ...
                100*sum(exclude)/nTraces );
        disp( to_row(find(exclude)) );
    end
    
    
    % Save mean parameters of prior distribution as new model parameters.
    % The EB prior is essentially the distribution that encompasses all of
    % the observed parameter values from all of the traces.
    result = selfanalysis(nStates);
    
    optModels = model;
    optModels.mu = result.prior.mu;

    u_a = 0.5 .* result.prior.nu;
    u_b = 0.5 ./ result.prior.W;
    optModels.sigma = (u_a./u_b) .^ -0.5;  %state_stdev.m
    
    A = bsxfun( @rdivide, result.prior.A, sum(result.prior.A,2) );
    optModels.rates = logm(A) / (dt/1000);
    if any(optModels.rates(:)<0)
        optModels.rates = (A-eye(size(A))) / (dt/1000);  %eq. 170
        disp('logm(A) return negative rates');
    end
    
    
    % Re-assign idealization state assignments to fall into the closest
    % ensemble class by FRET value (so all traces are consistent).
    idl = zeros(size(data));  %state sequences
    
    for i=1:nTraces,
        traceStates = result.posterior(i).mu;
        v = result.viterbi(i).state;
        
        for k=1:numel(traceStates),
            d = abs(optModels.mu-traceStates(k));
            v(v==k) = find(d==min(d),1,'first');
        end
        
        idl(i,1:numel(v)) = v;
    end
    
    % 
    idlTotal = zeros(origSize);
    idlTotal( ~params.exclude, :) = idl(:,1:origSize(2));
    
    LL = L(end);
    close(wbh);
end




function analysis = init_analysis2(xall, num_states)
% Rewrite of the ebFRET internal MainWindow.init_analysis method
% to avoid references to the GUI layer. Unlike the original, this is
% intended only to construct a new analysis object, rather than potentially
% update an old one.

nTraces = numel(xall);

analysis = struct('dim',{}, 'prior',{}, 'posterior',{}, 'expect',{}, 'viterbi',{});
e = repmat( struct('z',[], 'z1',[], 'zz',[], 'x',[], 'xx',[]), nTraces,1);
v = repmat( struct('state',[],'mean',[]), nTraces,1 );


for k = num_states
    % set number of states
    analysis(k).dim.states = k;
    
    % assign prior
    analysis(k).prior = ebfret.analysis.hmm.guess_prior(xall, k);
    
    % posterior defaults to prior
    analysis(k).posterior = repmat( analysis(k).prior, nTraces,1 );
    
    % Initialize other parameters
    analysis(k).viterbi = v;
    analysis(k).lowerbound = zeros(nTraces,1);
    analysis(k).expect = e;
end

end




function selfanalysis = run_vbayes(x, selfanalysis, nStates, varargin)
% Variational Bayes (vbFRET) algorithm for model optimization and idealization.
%    x is a cell array of FRET traces.
%    selfanalysis is the ebFRET internal structure for analysis state.
%    nStates is the model complexity to fit (may be a list)

    ip = inputParser();
    ip.StructExpand = true;
    ip.addParamValue('restarts', 0, @isscalar);       
    ip.addParamValue('threshold', 1e-6, @isscalar);       
    ip.addParamValue('max_iter', 100, @isscalar);       
    ip.addParamValue('soft_kmeans', false, @isscalar);       
    ip.parse(varargin{:});
    args = ip.Results;
    
    b_size = 24;  %traces per batch
    
    nTraces = numel(x);

    % loop over different analysis sets
    for a = nStates,
        
        % break up series into batches
        batches = {};
        b_start = 1:b_size:nTraces;
        for b = 1:(length(b_start)-1)
            batches{b} = b_start(b):(b_start(b+1)-1);
        end
        batches{end+1} = b_start(end):nTraces;
        
        % get prior
        u = selfanalysis(a).prior;
        
        % populate posterior field names if necessary
        if isempty(selfanalysis(a).posterior)
            selfanalysis(a).posterior = u([]);
        end
        
        % loop over time series
        for b = 1:length(batches)

            % environment for parfor loop
            ns = batches{b};

            [w(ns)] = selfanalysis(a).posterior(ns);
            
            % output variables for parfor loop
            posterior = {};
            expect = {};
            viterbi = {};
            lowerbound = {};
            restart = {};
            
            % parallel loop over batch of traces
            parfor n = ns
                % construct initial guesses for posterior parameters
                w0 = u([]);
                if ~isempty(w(n).mu)
                    % always add a restart based on last result if it exists
                    w0(end+1) = w(n);
                end
                if args.restarts > 0
                    % first restart is uninformative
                    w0(end+1) = ebfret.analysis.hmm.init_posterior(x{n}, u);
                end
                for r = 2:args.restarts
                    % next restarts are seeded with random draws from prior
                    w0(end+1) = ebfret.analysis.hmm.init_posterior(...
                                    x{n}, u, ...
                                    'draw_params', true, ...
                                    'soft_kmeans', args.soft_kmeans);
                end

                % run variational bayes for each restart
                vb = struct();
                for r = 1:length(w0)
                    [vb(r).w,vb(r).L,vb(r).E] = ...
                        ebfret.analysis.hmm.vbayes(x{n}, w0(r), u);
                end

                % determine best result
                L_max = vb(1).L(end);
                r_max = 1;
                for r = 2:length(vb)
                    if (vb(r).L(end) - L_max) > 1e-2 * args.threshold * abs(L_max) ...
                        || isnan(L_max)
                        r_max = r;
                        L_max = vb(r).L(end);
                    end
                end

                if ebfret.analysis.hmm.valid_prior(vb(r_max).w)
                    % keep best result
                    lowerbound{n} = vb(r_max).L(end);
                    posterior{n} = vb(r_max).w;
                    restart{n} = r_max + args.restarts - length(w0);

                    % calculate viterbi path
                    % (this might only be necessary on the last interation)
                    [viterbi{n}.state, viterbi{n}.mean] = ...
                        ebfret.analysis.hmm.viterbi_vb(vb(r_max).w, x{n});

                    % remap expected statistics
                    expect{n}.z = sum(vb(r_max).E.gamma(2:end,:), 1)';
                    expect{n}.z1 = vb(r_max).E.gamma(1, :)';
                    expect{n}.zz = squeeze(sum(vb(r_max).E.xi, 1));
                    expect{n}.x = vb(r_max).E.xmean(:);
                    expect{n}.xx = vb(r_max).E.xvar + vb(r_max).E.xmean.^2;
                else
                    % Invalid posterior. Failed?
                    warning('ebFRET:InvalidPosterior', ...
                            'Unable to obtain a valid VB estimate for time series %d. Resetting posterior.', n);
                    lowerbound{n} = 0;
                    posterior{n} = u;
                    expect{n} = struct('z', [], 'z1', [], 'zz', [], 'x', [], 'xx', []);
                    viterbi{n} = struct('state', [], 'mean', []);
                    restart{n} = -1;
                end
            end

            % store batch results
            ns = ns(~cellfun(@isempty, posterior(ns)));
            selfanalysis(a).posterior(ns) = [posterior{:}];
            selfanalysis(a).expect(ns) = [expect{:}];
            selfanalysis(a).viterbi(ns) = [viterbi{:}];
            selfanalysis(a).lowerbound(ns) = [lowerbound{:}];
            selfanalysis(a).restart(ns) = [restart{:}];

        end
    end
end




