%% Load the Data
%  load('')
close all

addpath('fig')

figure('Units','inch','Position',[0 0 6 6])

set(gcf,'color','w')

movegui(gcf,'center')


idx = (floor(chain.length*0.1):1:chain.length);
%% Donor Channel

axes1 = subplot(2,2,1);
p1    = plot(chain.i,(chain.mu_D(1,:)),...
             chain.i(idx),chain.mu_D(1,idx),...
             chain.i,(chain.mu_D(2,:)),...
             chain.i(idx),(chain.mu_D(2,idx)),'.-');
p1(2).Color = 'g';
p1(4).Color = 'g';      
l1    = line(axes1.XLim,chain.params.ground.mu_D(1).*[1 1],'linestyle','-','color','g','LineWidth',1);
line(axes1.XLim,chain.params.ground.mu_D(2).*[1 1],'linestyle','-','color','g','LineWidth',1);

leg = legend([p1(1) p1(2) l1(1)],'Burn-in','{State Level}','{True Escape Rate}');
leg.Location = 'northoutside'
leg.NumColumns = 3;
leg.Position=[0.320158108273749,0.954970499035015,0.33596837393969,0.020512820054323];
legend box off

ylabel('\mu^{D}_{\sigma_{1,2}}');



axes11 = subplot(2,2,2);
histogram(chain.mu_D(1,idx),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
hold on
histogram((chain.mu_D(2,idx)),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(4),'color'));
hold off
ylim(axes1.YLim);
ytemp = linspace(axes1.YLim(1),axes1.YLim(end),100);
line(gampdf(ytemp,chain.params.phi_D,(chain.params.psi_D/chain.params.phi_D)),ytemp,'linestyle','--','color','m','linewidth',2);
line(axes11.XLim,chain.params.ground.mu_D(1).*[1 1],'linestyle','-','color','g','LineWidth',1);
line(axes11.XLim,chain.params.ground.mu_D(2).*[1 1],'linestyle','-','color','g','LineWidth',1);

%% Acceptor Channel
axes1_A= subplot(2,2,3);
p1_A    = plot(chain.i,(chain.mu_A(1,:)),...
               chain.i(idx),chain.mu_A(1,idx),...
               chain.i,(chain.mu_A(2,:)),...
               chain.i(idx),(chain.mu_A(2,idx)),'.-');
p1_A(2).Color = 'r';
p1_A(4).Color = 'r';          
l1_A    = line(axes1_A.XLim,chain.params.ground.mu_A(1).*[1 1],'linestyle','-','color','r','LineWidth',1);
line(axes1_A.XLim,chain.params.ground.mu_A(2).*[1 1],'linestyle','-','color','r','LineWidth',1);

% legend([p1_A(1) p1_A(2) l1_A(1)],'Burn-in','{State Level}','{True Escape Rate}');
ylabel('\mu^{A}_{\sigma_{1,2}}');



axes11_A = subplot(2,2,4);
histogram(chain.mu_A(1,idx),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1_A(2),'color'));
hold on 
histogram((chain.mu_A(2,idx)),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1_A(4),'color'));
hold off
ylim(axes1_A.YLim);
ytemp = linspace(axes1_A.YLim(1),axes1_A.YLim(end),100);
line(gampdf(ytemp,chain.params.phi_A,(chain.params.psi_A/chain.params.phi_A)),ytemp,'linestyle','--','color','m','linewidth',2);
line(axes11_A.XLim,chain.params.ground.mu_A(1).*[1 1],'linestyle','-','color','c','LineWidth',1);
line(axes11_A.XLim,chain.params.ground.mu_A(2).*[1 1],'linestyle','-','color','c','LineWidth',1);

