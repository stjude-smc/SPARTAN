%% Load the Data
%  load('')
close all

addpath('fig')

figure('Units','inch','Position',[0 0 3 3])

set(gcf,'color','w')

movegui(gcf,'center')
s = size(IT);
chain.length = s(1);
chain.i      = (1:1:chain.length);
idx          = (floor(chain.length*0.1):1:chain.length);
%% Donor Channel

axes1 = subplot(4,2,1);
p1    = plot(chain.i,(mu_D(1,:)),...
             chain.i(idx),mu_D(1,idx),...
             chain.i,(mu_D(2,:)),...
             chain.i(idx),(mu_D(2,idx)),'.-');
l1    = line(axes1.XLim,params.ground.mu_D(1).*[1 1],'linestyle','-','color','g','LineWidth',1);
line(axes1.XLim,params.ground.mu_D(2).*[1 1],'linestyle','-','color','g','LineWidth',1);
ylabel('\mu^{D}_{\sigma_{1},sigma_{2}}');



axes11 = subplot(4,2,2);
histogram(mu_D(1,idx),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
hold on
histogram((mu_D(2,idx)),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
hold off
ylim(axes1.YLim);
ytemp = linspace(axes1.YLim(1),axes1.YLim(end),100);
line(gampdf(ytemp,params.phi_D,(params.psi_D/params.phi_D)),ytemp,'linestyle','--','color','m','linewidth',2);
line(axes11.XLim,params.ground.mu_D(1).*[1 1],'linestyle','-','color','g','LineWidth',1);
line(axes11.XLim,params.ground.mu_D(2).*[1 1],'linestyle','-','color','g','LineWidth',1);



%% Acceptor Channel
axes1_A= subplot(4,2,3);
p1_A    = plot(chain.i,(mu_A(1,:)),...
               chain.i(idx),mu_A(1,idx),...
               chain.i,(mu_A(2,:)),...
               chain.i(idx),(mu_A(2,idx)), '.-');
l1_A    = line(axes1_A.XLim,params.ground.mu_A(1).*[1 1],'linestyle','-','color','r','LineWidth',1);
line(axes1_A.XLim,params.ground.mu_A(2).*[1 1],'linestyle','-','color','r','LineWidth',1);
ylabel('\mu^{A}_{\sigma_{1},\sigma_{2}}');



axes11_A = subplot(4,2,4);
histogram(mu_A(1,idx),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1_A(2),'color'));
hold on 
histogram((mu_A(2,idx)),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1_A(2),'color'));
hold off
ylim(axes1_A.YLim);
ytemp = linspace(axes1_A.YLim(1),axes1_A.YLim(end),100);
line(gampdf(ytemp,params.phi_A,(params.psi_A/params.phi_A)),ytemp,'linestyle','--','color','m','linewidth',2);
line(axes11_A.XLim,params.ground.mu_A(1).*[1 1],'linestyle','-','color','c','LineWidth',1);
line(axes1_A.XLim,params.ground.mu_A(2).*[1 1],'linestyle','-','color','r','LineWidth',1);

%% Donor Channel BACKGROUND
axes1_B = subplot(4,2,5);
p1    = plot(chain.i,(mu_back_D(1,:)),...
    chain.i(idx),mu_back_D(1,idx),'.-');
l1    = line(axes1_B.XLim,params.ground.mu_back_D(1).*[1 1],'linestyle','-','color','g','LineWidth',1);
ylabel('\mu^{D}_{{back}}');



axes11_B = subplot(4,2,6);
histogram(mu_back_D(1,idx),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
ylim(axes1_B.YLim);
ytemp = linspace(axes1_B.YLim(1),axes1_B.YLim(end),100);
line(gampdf(ytemp,params.chi_D,(params.nu_D/params.chi_D)),ytemp,'linestyle','--','color','m','linewidth',2);
line(axes11_B.XLim,params.ground.mu_back_D(1).*[1 1],'linestyle','-','color','g','LineWidth',1);





%% Acceptor Channel BACKGROUND
axes1_BA= subplot(4,2,7);
p1_A    = plot(chain.i,(mu_back_A(1,:)),...
    chain.i(idx),mu_back_A(1,idx),'.-');
l1_A    = line(axes1_BA.XLim,params.ground.mu_back_A(1).*[1 1],'linestyle','-','color','r','LineWidth',1);
ylabel('\mu^{A}_{{back}}');



axes11_BA = subplot(4,2,8);
histogram(mu_back_A(1,idx),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1_A(2),'color'));
ylim(axes1_BA.YLim);
ytemp = linspace(axes1_BA.YLim(1),axes1_BA.YLim(end),100);
line(gampdf(ytemp,params.chi_A,(params.nu_A/params.chi_A)),ytemp,'linestyle','--','color','m','linewidth',2);
line(axes11_BA.XLim,params.ground.mu_back_A(1).*[1 1],'linestyle','-','color','c','LineWidth',1);


