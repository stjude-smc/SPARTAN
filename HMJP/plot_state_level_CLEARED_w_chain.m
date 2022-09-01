%% Load the Data
%  load('')
close all

addpath('fig')

figure(22)

figure('Units','inch','Position',[0 0 4 4])

set(gcf,'color','w')

movegui(gcf,'center')

s= size(chain.IT);
chain.length = s(1);
idx = (floor(0.1*chain.length):1:chain.length);
chain.i = 1:1:chain.length;

%% Donor Channel
axes1 = subplot(3,2,1);
p1    = plot(chain.i,(chain.mu_D(1,:)),...
            chain.i(idx),chain.mu_D(1,idx),...
            chain.i,(chain.mu_D(2,:)),...
            chain.i(idx),(chain.mu_D(2,idx)),...
            chain.i,(chain.mu_A(1,:)),...
            chain.i(idx),chain.mu_A(1,idx),...
            chain.i,(chain.mu_A(2,:)),...
            chain.i(idx),(chain.mu_A(2,idx)),'.-');
ylabel('Emission rate\newline(photons/s)');
p1(2).Color = 'g';
p1(4).Color = 'g';
p1(6).Color = 'r';
p1(8).Color = 'r';
p1(1).Color = 'b';
p1(3).Color = 'b';
p1(5).Color = 'b';
p1(7).Color = 'b';
l1    = line(axes1.XLim,chain.params.ground.mu_D(1).*[1 1],'linestyle','-','color','g','LineWidth',1);
l1    = line(axes1.XLim,chain.params.ground.mu_D(2).*[1 1],'linestyle','-','color','g','LineWidth',1);
l1    = line(axes1.XLim,chain.params.ground.mu_A(1).*[1 1],'linestyle','-','color','r','LineWidth',1);
l1    = line(axes1.XLim,chain.params.ground.mu_A(2).*[1 1],'linestyle','-','color','r','LineWidth',1);
axes11 = subplot(3,2,2);

ytemp = linspace(axes1.YLim(1),axes1.YLim(end),50);

histogram(chain.mu_D(1,idx),ytemp,'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
hold on
histogram((chain.mu_D(2,idx)),ytemp,'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(4),'color'));
histogram(chain.mu_A(1,idx),ytemp,'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(6),'color'));
histogram(chain.mu_A(2,idx),ytemp,'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(8),'color'));
ylim(axes1.YLim);
ytemp = linspace(axes1.YLim(1),axes1.YLim(end),100);
line(gampdf(ytemp,chain.params.phi_D,(chain.params.psi_D/chain.params.phi_D)),ytemp,'linestyle','--','color','m','linewidth',2);
l1    = line(axes11.XLim,chain.params.ground.mu_D(1).*[1 1],'linestyle','-','color','g','LineWidth',1);
l1    = line(axes11.XLim,chain.params.ground.mu_D(2).*[1 1],'linestyle','-','color','g','LineWidth',1);
l1    = line(axes11.XLim,chain.params.ground.mu_A(1).*[1 1],'linestyle','-','color','r','LineWidth',1);
l1    = line(axes11.XLim,chain.params.ground.mu_A(2).*[1 1],'linestyle','-','color','r','LineWidth',1);
%% Donor Channel BACKGROUND
axes1_B = subplot(3,2,3);
p1    = plot(chain.i,(chain.mu_back_D(1,:)),...
    chain.i(idx),chain.mu_back_D(1,idx),'.-');
p1(2).Color = 'g';
p1(1).Color = 'b';
l1    = line(axes1_B.XLim,chain.params.ground.mu_back_D(1).*[1 1],'linestyle','-','color','g','LineWidth',1);

% legend([p1(1) p1(2)],'Burn-in','{State Level}');
ylabel('Background\newlineemission rate\newline(photons/s)');



axes11_B = subplot(3,2,4);
ytemp = linspace(axes1_B.YLim(1),axes1_B.YLim(end),50);

histogram(chain.mu_back_D(1,idx),ytemp,'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
l1    = line(axes11_B.XLim,chain.params.ground.mu_back_D(1).*[1 1],'linestyle','-','color','g','LineWidth',1);
ylim(axes1_B.YLim);
line(gampdf(ytemp,chain.params.chi_D,(chain.params.nu_D/chain.params.chi_D)),ytemp,'linestyle','--','color','m','linewidth',2);


%% Acceptor Channel BACKGROUND
axes1_BA= subplot(3,2,5);
p1_A    = plot(chain.i,(chain.mu_back_A(1,:)),...
    chain.i(idx),chain.mu_back_A(1,idx),'.-');
p1_A(2).Color = 'r';
p1_A(1).Color = 'b';

% legend([p1_A(1) p1_A(2)],'Burn-in','{State Level}');
ylabel('Background\newlineemission rate\newline(photons/s)');
xlabel('MCMC iteration');
l1    = line(axes1_BA.XLim,chain.params.ground.mu_back_A(1).*[1 1],'linestyle','-','color','r','LineWidth',1);



axes11_BA = subplot(3,2,6);
ytemp = linspace(axes1_BA.YLim(1),axes1_BA.YLim(end),50);

histogram(chain.mu_back_A(1,idx),ytemp,'normalization','pdf','Orientation','horizontal','Facecolor',get(p1_A(2),'color'));
hold on
l1    = line(axes11_BA.XLim,chain.params.ground.mu_back_A(1).*[1 1],'linestyle','-','color','r','LineWidth',1);

ylim(axes1_BA.YLim);
line(gampdf(ytemp,chain.params.chi_A,(chain.params.nu_A/chain.params.chi_A)),ytemp,'linestyle','--','color','m','linewidth',2);
xlabel('Post. prob. distr.');


