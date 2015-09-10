function [energies] = elandscape(rates, pi, sampleLabels, stateLabels)
%ELANDSCAPE  Creates an array of contour, pop hist, and TD plots
% 
%   ELANDSCAPE(RATES, PI, sampleLabels, rateLabels)
%   Uses transition RATES from QuB (using calc_rates) and state occupancy
%   probabilities (from PercentTime) to construct an energy landscape plot.
%
%   RATES will have to be ordered in a specific way to make sense:
%   A->B, B->A, B->C, C-B, C->A, A-C  (forming a cycle)
%   Same order as in Munro et al 2007 (James' paper).

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% TODO: Use actual state probabilities for del-G standard to prevent
%   error accumulation.  G_dagger = AVG[ G_left G_right ];

%% INITIALIZE & PROCESS FUNCTION PARAMETERS

% Prase optional arguments
if nargin<3
    stateLabels  = {'C', 'H1', 'H2', 'C'};
end


% Reorder rates so that it follows the progression
% across the energy landscape:
%                 1      2      3      4      5       6
% INPUT from QuB: C->H1, H1->C, C->H2, H2->C, H1->H2, H2->H1
%
% OUTPUT for energy landscape: C->H1, H1->C, H1->H2, H2->H1, H2->C, C->H2
% energy landscape: C->H1->H2->C
% 
rates = rates(:,[1 2 5 6 4 3]);


%pi is in order of FRET value:
% H2, H1, C (in perecnt, 0-100)


% Expecting tRNA-tRNA FRET model (for now - generalize later)
[nSamples,nRates] = size(rates);

assert( size(rates,2)==6, 'unsupported model?' );


if nSamples ~= size(pi,1) || nRates/2 ~= size(pi,2),
    disp('DWT does not match rates');
end



%---- USER TUNABLE PARAMETERS ----

%---------------------------------




%% Calculate energies along coordinate  (v = h/(Kb*T))
% Source for constants = wikipedia
v = 6.62606896e-34 / (1.3806504e-23 * 298);  %sec
RT = 298 * 8.314472 / 1000;  %kJ/mol

% Calculate barrier energies from rates
Gdagger = -RT * log( rates*v );


% The old way using rates only
if ~exist('pi','var')
    % successively find peak and well energies by adding and
    % subtracting barrier energies
    Gdagger(:,2:2:end) = -Gdagger(:,2:2:end);
    energies = [zeros(nSamples,1) cumsum(Gdagger,2)];

% The "right way" to do it
else
    % Calculate stable state energies (relative to C)
    Gstandard = zeros(nSamples,4);
    Gstandard(:,2) = -RT*log( pi(:,2)./pi(:,3) );
    Gstandard(:,3) = -RT*log( pi(:,1)./pi(:,3) );

    energies = zeros(nSamples, 7);
    energies(:,1:2:end) = Gstandard;
    
    % Calculate transition state energies
    Gdagger2 = zeros(nSamples,3);
    Gdagger2(:,1) = Gdagger(:,1) + ( Gdagger(:,2)-Gstandard(:,2) );
    Gdagger2(:,2) = ( Gdagger(:,3)-Gstandard(:,2) ) + ( Gdagger(:,4)-Gstandard(:,3) );
    Gdagger2(:,3) = Gdagger(:,6) + ( Gdagger(:,5)-Gstandard(:,3) );
    
%     Gdagger2(:,1) = Gdagger(:,1);
%     Gdagger2(:,3) = Gdagger(:,6);
    
    energies(:,2:2:end) = Gdagger2/2;
end



if nargin<3
    return;
end


%% Plot the energy landscape

% Make cut in y-axis for better readability.
energies(:,2:2:end) = energies(:,2:2:end)-60;


% figure();
cla; hold on;

% First get colors for plot
j = colormap(jet);  nj = size(j,1);
colors = j( 1:floor(nj/nSamples):nj,: );


% solid points at calculated positions with the rest filled in with
% a simple spline to fill in the rest


% Plot dots for known points
x  = 1:(nRates+1);
xx = 1:0.01:x(end);
yy2 = zeros( nSamples, numel(xx) );

for i=1:nSamples,
    y = [0 energies(i,:) 0];
    yy2(i,:) = interp1( [0 x nRates+2], y, xx, 'pchip' );
end

% plot( xx, energies );
set(gca,'ColorOrder',colors(2:end,:));
plot( xx,yy2 );


% add circles for barriers and lines for wells
% plot( energies', 'o', 'MarkerSize',8 );
plot( [1,3,5,7], energies(:,1:2:end)', 'o', 'MarkerSize',8 );
plot( [2,4,6], energies(:,2:2:end)', 'o', 'MarkerSize',8 );


% Add labels etc...
if exist('sampleLabels','var')
    legend( sampleLabels );
end

%
set(gca, 'XTick', 1:7);

confLabel = cell( 1, nRates+1 );
confLabel(1:2:nRates+1) = stateLabels;
confLabel(2:2:nRates+1) = {'<->'};
set(gca, 'XTickLabel', confLabel );

xlabel('Conformation Coordinate');
ylabel( 'Energy (kJ/mol)' );
axis( [1 7 -1 12] );
grid on; box on;


% Create and plot spline fit







