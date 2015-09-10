function BWplot_dists(mus,sigmas,ps,data)
%
% Usage: BWplot_dists(mus,sigmas,ps,<data>)
%
% Makes nice histogram of model predictions versus data (if data is included)
%
% DAB 2008.3.30

%   Copyright 2008-2015 Cornell University All Rights Reserved.


RES = 0.01;

Nstates = length(mus);
[a b] = min(mus);
minval = mus(b)-sigmas(b)*4;
[a b] = max(mus);
maxval = mus(b)+sigmas(b)*4;

x = minval:RES:maxval;

dtot = zeros(1,length(x));
figure
if nargin > 3
  hist(data,minval:0.01:maxval)
  mult = length(data)*0.01;
else
  mult = 1;
end

hold on
cord = 'bcgrmbgcrm';
for n = 1:Nstates
  d = mult*ps(n)/sigmas(n)/sqrt(2*pi) * exp(-(x-mus(n)).^2/(2*sigmas(n)^2));
  dtot = dtot + d;
  plot(x,d,cord(n))
end

plot(x,dtot,'k')
axis([minval maxval 0 max(dtot)*1.1])
set(gca,'LineWidth',1.5,'FontSize',12)
