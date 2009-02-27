function [data,idl] =  simulate( Q, p0, kBleach, model, nTraces, traceLen, ...
                                 framerate, seed )
% SIMULATE   MCMC
%
%  
simFramerate = 1000; %1ms
binFactor = simFramerate/framerate;

tic;

mu = model(:,1);
sigma = model(:,2);
nStates = length(mu);

% Construct transition probability matrix (A) for 1ms simulation
A = Q./simFramerate;
A( find(eye(nStates)) ) = 1-sum(A,2);


data = zeros( nTraces, traceLen );  %FRET data at every time point
idl  = zeros( nTraces, traceLen*binFactor );  %state assignment at every time point


% Initialize the random number generator
if ~exist('seed','var')
    seed  = 1272729;
%     seed  = 1272720;
end

rand('twister',seed);
randn('state',seed);

% Generate a set of uniform random numbers for choosing states
rset = rand(nTraces,traceLen*binFactor);
  
h = waitbar(0,'Simulating...');
  
%--- Generate noiseless data at 1 ms
for i=1:nTraces,
    
    % Choose the initial state
    idl(i,1) = find( rset(i,1) <= cumsum(p0), 1, 'first' );
    
    % Choose each successive state    
    for j=2:traceLen*binFactor,
        choices = A( idl(i,j-1), : );
        idl(i,j) = find( rset(i,j) < cumsum(choices), 1, 'first' );
    end
    
    assert( all(idl(i,:)>0) );
    
    waitbar(i/(nTraces*2),h)
end



%--- Chop each trace to simulate exponential photobleaching
if kBleach > 0,
    pbTimes = -log(rand(1,nTraces))./kBleach;
    pbTimes = ceil(pbTimes*simFramerate);

    % remove times less than 15f*25ms= 375ms
    % simulating the filtering in autotrace

    for i=1:nTraces,
        idl(i,pbTimes(i):end) = 1;  % set the state to blinking
    end
end



%-- Add noise to the data
for i=1:nTraces,
    
    % Generate and and time-average the noiseless trace
    trace = mu( idl(i,:) );
    trace = binFretData( trace, binFactor );
    
    % Add gaussian noise
%     trace( idl(i,:)==1 ) = trace( idl(i,:)==1 ) + sigma(1)*randn(1,traceLen);
%     trace( idl(i,:)==2 ) = trace( idl(i,:)==2 ) + sigma(2)*randn(1,traceLen);
%     trace( idl(i,:)==3 ) = trace( idl(i,:)==3 ) + sigma(3)*randn(1,traceLen);
%     trace( idl(i,:)==4 ) = trace( idl(i,:)==4 ) + sigma(4)*randn(1,traceLen);
    trace = trace + sigma(end)*randn(1,traceLen);
    
    data(i,:) = trace;
    
    waitbar( (i+nTraces)/(2*nTraces), h)
end
close(h);

disp(toc);
end %FUNCTION simulate




function output = binFretData( trace, factor )

assert( factor >= 1 );
assert( floor(factor)==factor, 'binFretData: only integer values accepted' );

traceLen = numel(trace);

nFrames2 = floor(traceLen/factor);
output = zeros( 1,nFrames2 );

for j=1:nFrames2

    s = factor*(j-1) +1;
    e = factor*j;

    output( j ) = mean( trace(s:e) );

end


end %FUNCTION binFretData
