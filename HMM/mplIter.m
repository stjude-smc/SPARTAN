function [LL,dLL] = mplIter(data, dt, p0, classidx, params)
% Maximum Point Likelihood algorithm (MPL)
%
%   [LL,dLL] = mplIter(DATA, DT, RATES, p0, CLASSIDX, PARAMS)
%   Runs one iteration of the Maximum Point Likelihood algorithm (MPL).
%   DATA are experimental observations (one molecule trace per row).
%   DT is the experimental sampling interval of the data in seconds.
%   P0 is a vector of initial state probabilities.
%   CLASSIDX gives the class number of each state.
%
%   PARAMS is a vector of all optimizable parameters arranged as follows:
%     [mu1 mu2 ... stdev1 stdev2 ... k12 k13 .. k21 k23 .. k31 k32 .. ...].
%
%   LL is the log likelihood for current parameters: log[ P(data|model) ].
%   dLL is the gradient of the LL function w/r/t each element in PARAMS.
%   Both values are negated to facilitate use with minimizers (e.g., fmincon).
%  
%   See Qin et al (2000) Biophys J 79, pg. 1915-1927 for algorithm details.
%
%   See also: mplOptimize, bwOptimize, milOptimize, batchKinetics.
%

%   Copyright 2018 Cornell University All Rights Reserved.


%% Process input arguments
narginchk(5,5);
nargoutchk(0,2);

nStates = numel(p0);
p0 = reshape(p0, 1, nStates);
if isvector(data)
    data = reshape(data, 1, numel(data));
end
nTraces = size(data,1);

% Unpack parameters from fminunc input vector
mu     = params( 1:nStates );
sigma2 = params( nStates   + (1:nStates) );
rates  = params( 2*nStates + 1:end );
sigma = sqrt(sigma2);

Q = zeros(nStates);
I = logical(eye(nStates));
Q(~I) = rates;
Q(I) = -sum(Q,2);

% Calculate transition probability matrix (A) using spectral matrices.
% See pg. 615 of 1995 book "Single Channel Recording". See eq. 6 in Qin 2000.
% expm(Q*dt) would be simpler, but we need the spectral matrices anyway.
[righteig,eigval] = eig(Q);
eigval = diag(eigval);
transitionProb = zeros(nStates);
spectralMatrix = cell(nStates,1);

for s=1:nStates
    lefteig = righteig^-1;
    spectralMatrix{s} = righteig(:,s) * lefteig(s,:);
    transitionProb = transitionProb  +  exp(eigval(s)*dt) * spectralMatrix{s};
end


%% Calculate log likelihood and partial derivatives for each trace.
LL = 0;
dLL_mu     = zeros(1, numel(mu));
dLL_sigma2 = zeros(1, numel(mu));
dLL_A      = zeros(nStates);

for n=1:nTraces
    % Remove frames after donor photobleaching (which are precisely zero).
    trace = data(n,:);
    bleachFrame = find( trace==0, 1, 'first' );
    trace = trace(1:bleachFrame-1);
    nFrames = numel(trace);
    if nFrames<10, continue; end  %skip extremely short traces
    
    % Calculate emmission probabilities at each timepoint.
    % The factor of 100 ensures that observation probabilities are always less than
    % 1.0 (so log likelihood sign doesn't invert). This is the effective number of
    % bins for the distribution.
    observProb = zeros(nFrames, nStates);  %B matrix
    for i=1:nStates
        class = classidx(i);
        observProb(:,i) = normpdf( trace, mu(class), sigma(class) );
    end
    observProb = observProb+eps;  %avoid underflow to zero.


    %% Calculate forward/backward probabilities and log likelihood
    % Values are repeatedly scaled to prevent underflow of machine precision
    % NOTE: X * diag(Y) = X .* Y if X and Y are vectors of same size/orientation.
    % NOTE: this is the slowest step because of the direct iteration over frames.
    %       Consider implementing only this part as a mex function.
    % FIXME: identical to BWtransition. should combine.

    % Forward probabilities
    % alpha(t,i) = P( observations 1..t & state(t)=i | model )
    alpha = zeros(nFrames, nStates);
    nrm = zeros(nFrames,1);

    alpha(1,:) = p0 .* observProb(1,:);
    nrm(1) = 1./sum(alpha(1,:));
    alpha(1,:) = alpha(1,:) * nrm(1);

    for t=2:nFrames
        %alpha(t,j) = SUM_i( alpha(t-1,i)*A(i,j) ) * B(t,j)
        alpha(t,:) = (alpha(t-1,:) * transitionProb) .* observProb(t,:);

        nrm(t) = 1./sum(alpha(t,:));  %normalization prevents float underflow
        alpha(t,:) = alpha(t,:) * nrm(t);
    end

    LL = LL + -sum( log(nrm) );

    % Backward probabilities.
    % beta(t,i) = P( observations t+1..end | state(t)=i & model )
    beta = zeros(nFrames, nStates);
    beta(nFrames,:) = 1;

    for t=nFrames-1:-1:1
        %beta(t,i) = SUM_j(  A(i,j) * B(t+1,j) * beta(t+1,j)  )
        beta(t,:) = transitionProb * (observProb(t+1,:) .* beta(t+1,:))' * nrm(t);
    end

    % gamma(t,i) = P( state(t)=i | all obseravtions & model )
    % SUM_t(gamma) is the expected number of times each state is occupied.
    gamma = alpha .* beta;
    gamma = bsxfun( @rdivide, gamma, sum(gamma,2) );  %normalized at each time t


    %% Calculate partial derivatives of LL (see eq. 12-19)
    
    % Gradient of initial probabilities can be calculated like this, but in 
    % practice optimizing this parameter isn't worth the effort. Instead, we 
    % (plan to) calculate it from the final backward probabilities beta(:,1). FIXME
    % dLL_p0 = beta(1,:) .* observProb(1,:);   %eq. 12

    for i=1:nStates
        for j=1:nStates
            % Gradient of log likelihood by a specific rate constant (k_i,j)
            % This is not in the paper -- derived by me instead. Equivalent?
            %if i==j, continue; end
            %dA_dk = derivative_A_by_rate( i,j, spectralMatrix, eigval, dt);
            %for t=1:nFrames-1
            %    dLL_k(i,j) = dLL_k(i,j) + (alpha(t,:) * dA_dk * diag(observProb(t+1,:)) * beta(t+1,:)');
            %end

            % Gradient of transition probabilities from forward/backward variables.
            % sum_t[ alpha_t(j) * beta_t+1(j) * Bt+1(j) ]. See eq. 13
            dLL_A(i,j) = dLL_A(i,j) + sum(  alpha(1:end-1,i) .* beta(2:end,j) .* observProb(2:end,j)  );
        end

        class = classidx(i);
        dev = trace-mu(class);
        dLL_mu(i)     = dLL_mu(i)     + dev * gamma(:,i) ./ sigma2(class);
        dLL_sigma2(i) = dLL_sigma2(i) + ( -1/(2*sigma2(class)) + dev.^2/(2*sigma2(class)^2) ) * gamma(:,i);
        %dLL_sigma(i) = dLL_sigma(i)  + ( -1/sigma(class) + dev.^2/sigma(class)^3 )  *  gamma(:,i);
    end

end %for each trace


% Partial derivative of log likelihood by each rate constant using the chain
% rule. See eq. 19 (treating A as a vector to avoid tensor multiplication).
dLL_k = zeros(nStates);
for i=1:nStates
    for j=1:nStates
        if i==j, continue; end
        dA_dk = derivative_A_by_rate( i,j, spectralMatrix, eigval, dt);
        dLL_k(i,j) = dLL_A(:)' * dA_dk(:);
    end
end

% Combine partial derivatives for all parameters into a single vector for the
% optimizer. dLL_k includes derivatives only for off diagonals since the
% diagonal elements are implied.
% Constant factors comes form the normalization for normpdf above.
dLL = [dLL_mu dLL_sigma2 dLL_k(~I)'];

% dLL = dLL*100;

% Return opposite of LL and dLL since to convert maximizing LL into minimizing
% -LL for use with fminunc/fmincon.
LL = -LL;
dLL = -dLL;


end  %function mplIter



%%
function output = derivative_A_by_rate( idxIn,idxOut, spectralMatrix, eigval, dt)
% Derivative of the probability matrix A w/r/t a specific rate constant
% See eq. 17 of Qin 2000.
% NOTE: since only two delements of dQdk are nonzero (and one is the opposite of
% the other), it should be possible to simplify this calculation further.

nStates = numel(eigval);
output = 0;

% Derivative of rate matrix Q by a specific rate constant (k_i,j) is 1 at the
% corresponding location in Q, -1 on the diagonal, and zero otherwise.
dQdk = zeros(nStates);
dQdk(idxIn,idxOut) = 1;
dQdk(idxIn,idxIn) = -1;

% % Taylor series approximation
% for n=1:20
%     temp = 0;
%     for k=0:n-1
%         temp = temp + Q^k * dQdk * Q^(n-1-k);
%     end
%     output = output + ( dt^n/factorial(n) )  *  temp;
% end

% Spectral expansion (exact)
for i=1:nStates,
    for j=1:nStates
        if i==j
            temp = dt * exp(eigval(i)*dt);
        else
            temp = ( exp(eigval(i)*dt) - exp(eigval(j)*dt) )  /  ( eigval(i) - eigval(j) );
        end
        output = output + spectralMatrix{i} * dQdk * spectralMatrix{j} * temp;
    end
end

end  %function derivative_A_by_rate

