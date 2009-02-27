function B = constructB( obs, mus, sigmas )

NT = length(obs);
nStates = length(mus);

% Calculate emmission probabilities at each timepoint
% if the emmission came from each of the posssible states
C = 0.3989422804014327./sigmas; %gaussian leading coefficient
% C = 1./sqrt(2*pi*sigmas); %gaussian leading coefficient
D = (2*sigmas.^2);  %exponential denominator

B = zeros( NT,nStates );
for i=1:nStates
    B(:,i) = C(i) .* exp(-  ((obs-mus(i)).^2) ./ D(i)  );
end
