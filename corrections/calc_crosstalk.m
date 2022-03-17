function [output,n,x,f] = calc_crosstalk(input)
% CALC_CROSSTALK  estimate fraction of donor to acceptor spectral overlap.
%
%   CT = calc_crosstalk( INPUT ) estimates the fraction of donor fluorescence
%   appearing on the acceptor channel (CT; spectral crosstalk).
%   INPUT may be a TracesFret object or the path to a .traces file.
%
%   [CT,N,X,F] = calc_crosstalk( INPUT ) also provides histogram data used
%   for finding ensemble average value (X=bin centers; N=histogram counts)
%   as well as the fit object (F).
%
%   ALGORITHM:
%   For each trace, average 15 frames just before donor bleaching and
%   divide the mean acceptor by mean donor value to get a crosstalk
%   estimate. Accumulate these into a histogram and fit to a Gaussian to
%   get ensemble-averaged crosstalk value.
%
%   CAVEATS:
%   The above approach assumes that the acceptor tends to photobleach
%   faster than the donor and/or there are many donor-only molecules.
%   If this is not the case, the algorithm may give unexpected results.
%
%   See also: crosstalkcorrect, scaleacceptor, gammacorrect, correctTraces.

%   Copyright 2014-2022 Cornell University All Rights Reserved.


%% Parameter values

% Histogram bin centers for ensemble crosstalk distribution.
x = -0.1:0.01:0.25;



%% Process input arguments
[output,n,f] = deal([]);

narginchk(0,1);
nargoutchk(0,4);

if nargin<1
    input = getFile('*.traces;*.rawtraces');
end
if isempty(input), return; end
if ischar(input)
    data = loadTraces(input);
elseif isa(input,'Traces')
    data = input;
else
    error('Invalid input');
end

% Verify input is simple 2-channel FRET
if isa(data,'TracesFret4')
    warning('Multi-color data is not supported');
elseif ~all( strcmpi(data.channelNames,{'donor','acceptor','fret'}) )
    error('Single-channel data is not supported');
end


%% Calculate crosstalk values
nTraces = size(data.donor,1);

% Calculate crosstalk values as the fraction of intensity on the acceptor
% channel relative to the donor, after acceptor photobleaching. If there is no
% crosstalk, the acceptor intensity should be zero (and crosstalk is also zero).
crosstalk = NaN( nTraces,1 ); %NaN corresponds to no data.

for i = 1:nTraces
    fret = data.fret(i,:);
    
    % Find the point at which the donor photobleaches.
    donor_lt =  find( fret~=0, 1,'last' );
    if isempty(donor_lt) || donor_lt<20,  continue;  end
    
    idx = donor_lt-16:donor_lt-2;  %indexes of frames to use for calculation.
    if any(fret(idx)==0),  continue;  end  %avoid traces with donor blinking
    
    c = mean( data.acceptor(i,idx) ) / mean( data.donor(i,idx) );
    if c>x(1)&&c<x(end), crosstalk(i)=c; end
end

% Remove NaN values (from traces where crosstalk could not be estimated).
fprintf( [mfilename ': using %d of %d traces (%.0f%%).\n\n'], sum(~isnan(crosstalk)),  ...
          numel(crosstalk), 100*sum(~isnan(crosstalk))/numel(crosstalk) );
crosstalk = crosstalk( ~isnan(crosstalk) );
assert( ~isempty(crosstalk), 'No usable traces for calculating crosstalk!' );

% Create a histogram of crosstalk values and fit to a Gaussian to get a
% precise ensemble-average value.
% TO DO: ignore very low occupancy state.
n = hist( crosstalk, x )/numel(crosstalk);
try
    f = fit( x',n', 'gauss2' );
    if abs(f.b1-f.b2)<0.02 || f.a1<0.01 || f.a2<0.01
        f = fit( x',n', 'gauss1' );
        output = f.b1;
    else
        output = min(f.b1,f.b2);
    end
    %figure; bar(x,n); hold on; plot(f);
catch
    output = mean(x);
    warning('Failed to fit crosstalk distribution; using mean instead');
end


end %function calc_crosstalk


