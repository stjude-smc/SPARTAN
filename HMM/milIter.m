function LL = milIter(dwt, dt, model, rates)
% milIter  Likelihood of dwell times given a model (MIL optimizer)
%
%   LL = milIter(DWT, DT, MODEL, RATES) calculates the log likelihood of
%   dwell-times given current model parameter values. DWT is a cell array
%   (one per trace) of dwell times (first column) and observed class
%   (second column). DT is the sampling interval in seconds. MODEL is a
%   QubModel object with the starting parameter values. RATES is a list of
%   rates constants passed from the optimizer, which does not include any 
%   rate constants that are fixed in the model.
%
%   The analytical form of the likelihood function is from Qin et al (1997)
%   PRSL/B 264, p. 375, which describes the MIL algorithm.
%
%   See also: dtAdjustedQ, milOptimize, batchKinetics.


narginchk(4,4);
nargoutchk(0,1);



%% Construct normalized rate matrix (Q) from input rates
% FIXME: this fixes all rates set to zero. They should only be
% automatically fixed like this if BOTH directions are zero, indicating
% there is no connection between the two states.
% FIXME: rateMask should be a function in QubModel.

Q = model.rates;
I = logical(eye(model.nStates));
rateMask = ~I & Q~=0 & ~model.fixRates;
Q(rateMask) = rates;

% Normalize so rows sum to 0
Q(I) = 0;
Q(I) = -sum(Q,2);

% Adjust rates to account for missed events (eQ in literature).
% NOTE: this seems to add significantly to the number of iterations needed to
% converge. Maybe something is wrong with this function?
Q = dtAdjustedQ(Q, 0.7*dt, model.class);

% Use equilibrium probabilities if initial probabilities not specified.
% See eq. 17 on pg. 597 (chapter 20) of "Single Channel Recording" (1995).
% U = ones(1, model.nStates);
% S = [ Q U' ];
% p0 = U * (S * S')^-1;



%% Construct spectral matrices for each submatrix of Q that describes rate
% constants for transitions within that class (Qaa for states w/i class a).
% Used for the calculation of exp(Qaa*t) and its partial derivatives.
% See eq. 89, pg. 616 of "Single Channel Recording" (1995).
eigvals = cell(model.nClasses, 1);
spectra = cell(model.nClasses, 1);

for a=1:model.nClasses
    statesInClass = model.class==a;
    [X,lambda] = eig( Q(statesInClass,statesInClass) );
    eigvals{a} = diag(lambda);
    
    Y = X^-1;
    for i=1:sum(statesInClass)
        spectra{a}{i} = X(:,i) * Y(i,:);  
    end
end



%% Calculate dwell observation probabilities for each dwell (see eq. 4)
% G_abt(i,j) = P( stay in state i of class a for t seconds and transition to state j of class b | current state is i).
LL = 0;

for traceID=1:numel(dwt)
    dwellClass = dwt{traceID}(:,1);     %class of each dwell ('a' in paper notation)
    dwellTimes = dwt{traceID}(:,2)*dt;  %observed time in each dwell in seconds
    nDwells = numel(dwellClass);
    
    anorm = ones(nDwells,1);  %normalization constants at each time
    alpha_k = to_row( model.p0(model.class==dwellClass(1)) );  %initial probabilities
    
    for k=1:nDwells
        aa = dwellClass(k);
        a = model.class==aa;  %states in current observed class
        
        % If this is the final dwell in the dark state, don't include it
        % in the analysis, since it adds no useful information.
        % FIXME: assumes dark state is the first class.
        if k==nDwells && aa==1
            anorm(k) = [];
            break;
        end
        
        % Calculate expm(Qaa*t)
        if sum(a)==1
           % Scalar version (no degenerate states)
           obsProb = exp( Q(a,a) * dwellTimes(k) );
        else
            % Spectral expansion of matrix exponential (see eq. 19)
            obsProb = 0;
            for i=1:numel( spectra{aa} )
                obsProb = obsProb + spectra{aa}{i} * exp( eigvals{aa}(i) * dwellTimes(k) );
            end
        end
        
        if k<nDwells
            b = model.class==dwellClass(k+1);
            obsProb = obsProb * Q(a,b);
        else
            % Termination:: marginalize over states from all other classes
            % since the target class of the final transition is unknown.
            obsProb = obsProb * Q(a,~a);
            obsProb = obsProb * ones( size(obsProb,2), 1 );
        end
        obsProb(obsProb==0)=eps;  %avoid underflow errors
        
        % Calculate normalized forward probabilities for log-likelihood
        alpha_k = alpha_k * obsProb;
        anorm(k) = 1/sum(alpha_k);
        alpha_k = alpha_k * anorm(k);

    end %for each dwell

    % Calculate -log likelihood from normalization coefficients.
    % The negative factor will make the minimizer find the maximum LL.
    LL = LL + sum( log(anorm) );
    
end %for each trace




end









