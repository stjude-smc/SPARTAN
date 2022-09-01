%% Load the Data
%  load('')
close all

addpath('fig')

figure('Units','inch','Position',[0 0 5 5])

set(gcf,'color','w')

movegui(gcf,'center')



m_length     = size(chain.Q);
chain.length = m_length(3);


idx = (floor(chain.length*0.2):1:chain.length);
chain.i = (1:1:chain.length);

axes1 = subplot(3,2,1);
p1    = plot(chain.i,-squeeze(chain.Q(1,1,:)),...
    chain.i(idx),-squeeze(chain.Q(1,1,idx)),'.-');
if isfield(chain.params,'ground')
l1    = line(axes1.XLim,chain.params.ground.escrate(1).*[1 1],'linestyle','-','color','c','LineWidth',1);
end
legg = legend([p1(1) p1(2) l1(1)],'Burn-in','{Escape Rate}','{True Escape Rate}','Orientation','horizontal');
legg.Location = 'northoutside';
legg.Position = [0.093157886748168,0.957426114569578,0.717368429041305,0.035789474233696];
legg.Box = 'off';
ylabel('q_{1} (1/s)');



axes11 = subplot(3,2,2);
histogram(-squeeze(chain.Q(1,1,idx)),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
ylim(axes1.YLim);
ytemp = linspace(axes1.YLim(1),axes1.YLim(end),100);
line(gampdf(ytemp,chain.params.eta/2,2/(chain.params.eta*chain.params.beta)),ytemp,'linestyle','--','color','m','linewidth',2);
if isfield(chain.params,'ground')
l1    = line(axes11.XLim,chain.params.ground.escrate(1).*[1 1],'linestyle','-','color','c','LineWidth',1);
end

axes2= subplot(3,2,3);
p1 = plot(chain.i,-squeeze(chain.Q(2,2,:)),...
    chain.i(idx),-squeeze(chain.Q(2,2,idx)),'.-');
    line(axes2.XLim,chain.params.ground.escrate(2).*[1 1],'linestyle','-','color','c','LineWidth',1);
 ylabel('q_{2} (1/s)');


axes21 = subplot(3,2,4);
histogram(-squeeze(chain.Q(2,2,idx)),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
ylim(axes2.YLim);
ytemp = linspace(axes2.YLim(1),axes2.YLim(end),100);
line(gampdf(ytemp,chain.params.eta/2,2/(chain.params.eta*chain.params.beta)),ytemp,'linestyle','--','color','m','linewidth',2);
if isfield(chain.params,'ground')
line(axes21.XLim,chain.params.ground.escrate(2).*[1 1],'linestyle','-','color','c','LineWidth',1);
end


axes3 = subplot(3,2,5);
p1 = plot(chain.i,-squeeze(chain.Q(2,2,:)),...
    chain.i(idx),-squeeze(chain.Q(2,2,idx)),'.-');
    line(axes3.XLim,chain.params.ground.escrate(3).*[1 1],'linestyle','-','color','c','LineWidth',1);
if isfield(chain.params,'ground')
line(axes3.XLim,chain.params.ground.escrate(3).*[1 1],'linestyle','-','color','c','LineWidth',1);
end
    ylabel('q_{3} (1/s)');

xlabel('MCMC Chain')

axes31 = subplot(3,2,6);
histogram(-squeeze(chain.Q(2,2,idx)),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
ylim(axes3.YLim);
ytemp = linspace(axes3.YLim(1),axes3.YLim(end),100);
line(gampdf(ytemp,chain.params.eta/2,2/(chain.params.eta*chain.params.beta)),ytemp,'linestyle','--','color','m','linewidth',2);
if isfield(chain.params,'ground')
line(axes31.XLim,chain.params.ground.escrate(3).*[1 1],'linestyle','-','color','c','LineWidth',1);
end
xlabel('PDF');
% 
