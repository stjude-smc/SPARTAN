function kl=klocking(rates)
%

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% Initialize rate values
kCH1  = rates(1);
kH1C  = rates(2);
kH1H2 = rates(3);
kH2H1 = rates(4);
kCH2  = rates(5);
kH2C  = rates(6);

kCB=0.05;
kCPB=0.05;

kH1B=0.05;
kH1PB=0.05;

kH2B=0.05;
kH2PB=0.05;

kBC=10;
kBH1=10;
kBH2=10;
kBPB=0;

kPBC=0.000001;
kPBH1=0.000001;
kPBH2=0.000001;
kPBB=0;

% Perform the calculation
phi=[kCH1/(kCH1+kCH2) kCH2/(kCH1+kCH2)];

kCC=-kCH1-kCH2-kCPB-kCB;
kH1H1=-kH1C-kH1H2-kH1B-kH1PB;
kH2H2=-kH2C-kH2H1-kH2B-kH2PB;
kBB=-kBC-kBH1-kBH2-kBPB;
kPBPB=-kPBC-kPBH1-kPBH2-kPBB;

K=[kH1H1 kH1H2 kH1C kH1B kH1PB;
   kH2H1 kH2H2 kH2C kH2B kH2PB;
   kCH1  kCH2  kCC  kCB  kCPB;
   kBH1  kBH2  kBC  kBB  kBPB;
   kPBH1 kPBH2 kPBC kPBB kPBPB];

[M,l]=eig(K(1:2,1:2));

lambda1=-l(1,1);
lambda2=-l(2,2);

N=M^-1;

A1=M(:,1)*N(1,:);
A2=M(:,2)*N(2,:);

N=(-K(1:2,1:2))^-1*K(1:2,3);

w1=phi*A1*K(1:2,3);
w2=phi*A2*K(1:2,3);

kl=(lambda2^2*lambda1^2)/(w1*lambda2^2+w2*lambda1^2);
