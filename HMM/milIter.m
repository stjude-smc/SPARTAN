function [LL,dLL] = milIter(dwt, dt, p0, classidx, rateMask, rates)
% Function is nested so we have access to the eigenvalues/spectral matrices.
% Unless specified, equations are from Qin et al (1997) PRSL/B 264 p. 375.
% 
%   Q is the matrix of rate constants normalized so that all rows sum to zero.
%   dwt is an Nx2 matrix of idealized class numbers and times in the two
%     columns, respectively.
%   classidx is an Nx1 vector of the class number associated with each state.
%

LL = 0;
dLL = zeros(size(rates));

% disp(rates);

nClasses = max(classidx);
nStates  = numel(classidx);
nTraces  = numel(dwt);



%% Construct normalized rate matrix (Q) from input rates
% Eq. 29 suggests a more correct way to do this, but I don't understand it.
Q = zeros(nStates);
Q(rateMask) = rates;
Q = Q+eps; %
Q( logical(eye(nStates)) ) = -sum(Q,2);  %normalize so rows sum to zero

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
% NOTE: expm() would work, but need spectral matrices anyway for derivatives.
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

for traceID=1:nTraces
    dwellClass = dwt{traceID}(:,1);     %class of each dwell ('a' in paper notation)
    dwellTimes = dwt{traceID}(:,2)*dt;  %observed time in each dwell in seconds
    nDwells = numel(dwellClass);
    obsProb = cell(nDwells, 1);
    
    for t=1:nDwells
        curclass = dwellClass(t);

        % Calculate expm(Qaa*t) from spectral expansion of Qaa (see eq. 19)
        %   = SUM_(state i in class a)[ A_a,i * exp(lambda_a,i*t) ]
        %Gt = exp( Q(curclass,curclass) * dwellTimes(t) ) + eps;  %non-degenerate equivalent
        Gt = eps;  %avoid zero probabilities with disallowed transitions.
        for i=1:numel( spectra{curclass} )
            Gt = Gt + spectra{curclass}{i} * exp( eigvals{curclass}(i) * dwellTimes(t) );
        end
        
        % Calculate G_ab(t) = expm(Qaa*t) * Qab
        if t<nDwells
            obsProb{t} = Gt * Q( classidx==curclass, classidx==dwellClass(t+1) );
        else
            % Termination case for final dwell: marginalize over all possible end 
            % states from all classes since the final transition is unknown.
            notcur = classidx~=curclass;
            obsProb{t} = Gt * Q(classidx==curclass, notcur) * ones(sum(notcur),1);
        end
    end %for each dwell
    
    
    % Forward probabilities and normalization
    % Elements of alpha and beta are all column vectors.
    % NOTE: indexing is the same as the paper, we just combine alpha 0 and 1.
    % FIXME: this will cause problems for calculating gradients.
    % NOTE: each alpha is a row vector of varying length.
    alpha = cell(nDwells,1);
    anorm = zeros(nDwells,1);  %normalization constants at each time

    alpha{1} = to_row( p0(classidx==dwellClass(1)) ) * obsProb{1};
    anorm(1) = 1./sum(alpha{1});
    alpha{1} = alpha{1}*anorm(1);

    for t=2:nDwells
        alpha{t} = alpha{t-1} * obsProb{t};  %eq. 10; 1xN * NxN = 1xN
        anorm(t) = 1./sum(alpha{t});
        alpha{t} = alpha{t} * anorm(t);
    end
    
    assert( ~any(isnan(log(anorm))) );

    % Calculate -log likelihood from normalization coefficients.
    % The negative factor will make the minimizer find the maximum LL.
    % NOTE: this is in units of log( probability per frame ) ??
    LL = LL + sum( log(anorm) );
    
end %for each trace



if nargout<2, return; end


% %% Backward probabilities. See eq. 12
% beta = cell(nDwells+1,1);
% beta{end} = ones( size(obsProb{t},2), 1 );
% 
% for t=nDwells:-1:1
%     beta{t} = obsProb{t} * beta{t+1} * anorm(t);
% end
% 
% % Calculate derivatives of likelihood function
% delQ = zeros(size(Q));  %derivatives of Q
% for a=1:nClasses
%     aidx = classidx==a;  %indexes of states in class a
%     
%     for b=1:nClasses
%         bidx = classidx==b;  %indexes of states in class b
%         
%         if a~=b
%             delQ(aidx,bidx) = calcDelQab(a,b);
% %         else
% %             delQ(aidx,bidx) = calcDelQaa(a);
%         end
%     end 
% end
% 
% dLL = delQ(rateMask);  %FIXME: only calculate these elements
% 

% 
%     %% Function to calculate the partial derivative of rate matrix Q
%     function dQab = calcDelQab(a,b)
%         dQab = zeros( sum(aidx), sum(bidx) );
% 
%         % Get indexes of dwells in class a that transition to class b.
%         dwellab = find( dwt(:,1:end-1)==a & dwt(:,2:end)==b );
%         lambdaa = lambda{a};
% 
%         % Outer sum over spectral matrix components
%         for ii=1:numel(lambdaa)
%             % Inner sum over dwells transitioning from class a to b
%             isum = 0;
%             for k=to_row(dwellab)
%                 dwelltime = dwt(k,2);
%                 isum = isum + alpha{k} * beta{k}' * exp(lambdaa(ii)*dwelltime);
%             end
%             dQab = dQab + A{a}{ii}' * isum;  %unclear if should be transposed!
%         end
% 
%     end %function delQab
% 
% 
%     function dQaa = calcDelQaa(a)
%         Na = sum(aidx);
%         dQaa = zeros( sum(aidx) );
% 
%         % Get indexes of dwells in class a that transition to class b.
%         dwella = find( dwt==a );
%         lambdaa = lambda{a};
% 
%         for ii=1:Na
%             for jj=1:Na
%                 dwellsum = 0;
%                 
%                 % For each dwell transitioning from class a to b...
%                 for k=to_row(dwella)
%                     dwelltime = dwt(k,2);
%                     
%                     if ii==jj
%                         gamma = dwelltime * exp( lambdaa(i) * dwelltime );
%                     else
%                         gamma = ( exp(lambdaa(ii)*dwelltime) - exp(lambdaa(jj)*dwelltime) )  ...
%                               /  (lambdaa(ii)-lambdaa(jj));
%                     end
%                     
%                     % From Qin 1996 -- different from Qin 1997???
%                     %Qtrans = Q( dwt(k,1), dwt(k+1,1) );
%                     dwellsum = dwellsum + alpha{k} * beta{k}' * gamma;
%                 end
%                 
%                 dQaa = dQaa + ( A{a}{ii}' * dwellsum * A{a}{jj}' );
%             end
%         end
% 
%     end %function delQab


    end









