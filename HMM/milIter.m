function LL = milIter(dwt, dt, p0, classidx, rateMask, rates)
% Function is nested so we have access to the eigenvalues/spectral matrices.
% Unless specified, equations are from Qin et al (1997) PRSL/B 264 p. 375.
% 
%   Q is the matrix of rate constants normalized so that all rows sum to zero.
%   dwt is an Nx2 matrix of idealized class numbers and times in the two
%     columns, respectively.
%   classidx is an Nx1 vector of the class number associated with each state.


narginchk(6,6);
nargoutchk(0,1);



%% Construct normalized rate matrix (Q) from input rates
nStates = numel(classidx);
Q = zeros(nStates);
Q(rateMask) = rates;
Q = Q+eps; %
Q( logical(eye(nStates)) ) = -sum(Q,2);  %normalize so rows sum to zero
% disp(Q);

% Adjust rates to account for missed events (referred to as eQ in literature).
% if ~isempty(dt)
%     Q = dtAdjustedQ(Q, 0.8*dt, classidx);
% end

% Use equilibrium probabilities if initial probabilities not specified.
% See eq. 17 on pg. 597 (chapter 20) of "Single Channel Recording" (1995).
% if isempty(p0)
%     U = ones(1, nStates);
%     S = [ Q U' ];
%     p0 = U * (S * S')^-1;
% end


%% Construct spectral matrices for each submatrix of Q that describes rate
% constants for transitions within that class (Qaa for states w/i class a).
% Used for the calculation of exp(Qaa*t) and its partial derivatives.
% See eq. 89, pg. 616 of "Single Channel Recording" (1995).
nClasses = max(classidx);

eigvals = cell(nClasses, 1);
spectra = cell(nClasses, 1);

for a=1:nClasses
    statesInClass = classidx==a;
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
    
    anorm = zeros(nDwells,1);  %normalization constants at each time
    alpha_k = to_row( p0(classidx==dwellClass(1)) );  %initial probabilities
    
    for k=1:nDwells
        a = classidx==dwellClass(k);
        
        % Calculate expm(Qaa*t)
        if sum(a)==1
           % Scalar version (no degenerate states)
           obsProb = exp( Q(a,a) * dwellTimes(k) ) + eps;
        else
            % Spectral expansion of matrix exponential (see eq. 19)
            obsProb = eps;  %avoid zero probabilities with disallowed transitions.
            for i=1:numel( spectra{a} )
                obsProb = obsProb + spectra{a}{i} * exp( eigvals{a}(i) * dwellTimes(k) );
            end
        end
        
        if k<nDwells
            b = classidx==dwellClass(k+1);
            obsProb = obsProb * Q(a,b);
        else
            % Termination:: marginalize over states from all other classes
            % since the target class of the final transition is unknown.
            obsProb = obsProb * Q(a,~a);
            obsProb = obsProb * ones( size(obsProb,1), 1 );
        end
        
        % Calculate normalized forward probabilities for log-likelihood
        alpha_k = alpha_k * obsProb;
        anorm(k) = 1/sum(alpha_k);
        alpha_k = alpha_k * anorm(k);

    end %for each dwell

    % Calculate -log likelihood from normalization coefficients.
    % The negative factor will make the minimizer find the maximum LL.
    % NOTE: this is in units of log( probability per frame ) ??
    assert( ~any(isnan(log(anorm))) );
    LL = LL + sum( log(anorm) );
    
end %for each trace



end









