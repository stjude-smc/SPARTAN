function [rates,errors,labels] = trnaRates(arg)
%CALC_RATES  Combines plots from multiple datasets
% 
%   Generates rates relevant to tRNA-tRNA FRET experiments using
%   rate data from QuB MIL together (4-state fully connected model)

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% TODO: ....


% No arguments given, ask user for a rates file
if nargin<1,
    [file path]=uigetfile('*.txt','Choose a file:');
    filename = strcat(path,file);
    qub_rates = load(filename);

% A filename was given, load it as a rate matrix
elseif isa(arg,'char')
    filename = arg;
    qub_rates = load(filename);

% A matrix was given, load it directly as the rate matrix
elseif isa(arg,'double')
    qub_rates = arg;
    
end

[nSamples,nRates] = size(qub_rates);



% This is the standard rate order established in Munro '07.
labels = { 'C->H1', 'H1->C', 'H1->H2', 'H2->H1', 'C->H2', 'H2->C', ...
           'k_locking', 'k_unlocking', 'kAP', 'kAA' };
disp('C->H1, H1->C, H1->H2, H2->H1, C->H2, H2->C, k_locking, k_unlocking, kAP, kAA');
% Assumes fully-connected 4 state model (PB,H2,H1,C=1,2,3,4)
% 1->2, 1->3, 1->4, 2->1, 2->3, 2->4, 
% 3->1, 3->2, 3->4, 4->1, 4->2, 4->3




% Translate into rates of physical model:
if nRates==6,
    rates = qub_rates;

elseif nRates==12,      % Raw rate list from qub or HMM scripts
    rates = qub_rates(:, [12,9, 8,5, 11,6 ] );
    errors = []; 
    
elseif nRates==24,  % Same as above, but with symmetric errors included
    rates  = qub_rates(:, 1:2:end);
    rates  = rates(:,  [12,9, 8,5, 11,6] );
    
    errors = qub_rates(:, 2:2:end);
    errors = errors(:, [12,9, 8,5, 11,6] );

else
    error('Invalid rate matrix (%d)', nRates);
    
end


% Calculate derived rates from the above data (collective tRNA motion)
special_rates = zeros( nSamples, 4 );

for i = 1:nSamples,
    special_rates(i,1) = klocking( rates(i,:) );
    special_rates(i,2) = rates(i,1) + rates(i,5);
    special_rates(i,3) = kAP( rates(i,:) );
    special_rates(i,4) = rates(i,2) + rates(i,3);
end

rates = [rates special_rates];

