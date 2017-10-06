function varargout = tirfProfile(input,bins)
% tirfProfile  Line scan intensity profiles from fluorescence intensity
%
%    [PROFILE,N,BINS] = TIRFPROFILE(FNAME) loads the .traces file in FNAME and
%    estimates the mean fluorescence intensity profile across the field of view.
%    PROFILE has distributions for both X and Y axis (first a second columns).
%    N is the number of molecules contributing to each bin of these histograms.
%    BINS are the bin centers. For best results, combine .rawtraces files from
%    several movies to use a large number of molecules.
%
%    TIRFPROFILE(FNAME) with no output arguments displays the result instead.
%
%    See also: gettraces.

%    Copyright 2017 Cornell University All Rights Reserved.



%% Process input arguments
narginchk(0,1);
nargoutchk(0,3);

if nargin<1
    input = getFile;
    if isempty(input), return; end  %user hit cancel
end
if ischar(input),
    data = loadTraces(input);
elseif isa(input,'Traces')
    data = input;
else
    error('Invalid input. Must be a filename or Traces object');
end



%% Calculate intensity profiles
x = [data.traceMetadata.donor_x];
y = [data.traceMetadata.donor_y];

if nargin<2,
    bins = round(linspace(1,max(x)-min(x),25));
end

% stats = traceStat(data);
% t = [stats.t];
t = mean(data.total(:,1:10),2); %much faster

% Sum intensities across X and Y
profile = zeros( numel(bins)-1, 2 );
N = zeros( numel(bins)-1, 2 );

for i=1:numel(bins)-1
    sel = x>=bins(i) & x<bins(i+1);
    profile(i,1) =  sum( t(sel) )/sum(sel);  %average intensity in this region
    N(i,1) = sum(sel);
    
    sel = y>=bins(i) & y<bins(i+1);
    profile(i,2) =  sum( t(sel) )/sum(sel);
    N(i,2) = sum(sel);
end

binCenters = (bins(1:end-1)+bins(2:end))/2;
output = {profile,N,binCenters};
[varargout{1:nargout}] = output{1:nargout};

if nargout>0, return; end



%% Plot intensity distributions
figure;
subplot(2,2,1);
plot( binCenters, profile(:,1) );
title('X profile');
ylim([0 1.05*max(profile(:))])

subplot(2,2,2);
plot( binCenters, profile(:,2) );
title('Y profile');
ylim([0 1.05*max(profile(:))])

subplot(2,2,3);
plot( binCenters, N(:,1) );
title('Number of Molecules');
ylim([0 1.05*max(N(:))])

subplot(2,2,4);
plot( binCenters, N(:,2) );
title('Number of Molecules');
ylim([0 1.05*max(N(:))])


end



