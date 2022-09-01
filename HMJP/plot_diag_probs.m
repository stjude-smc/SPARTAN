%% Load the Data
%  load('')

figure('Units','inch','Position',[0 0 3 3])

set(gcf,'color','w')

movegui(gcf,'center')


idx = (1:1:chain.length);

pi_mjp        = nan(chain.params.M,chain.params.M,(chain.length));

for j=1:(chain.length)   
    
pi_mjp(:,:,j) = (expm((chain.params.t_left(2)-chain.params.t_left(1))*squeeze(chain.Q(:,:,j)))) ;

    if isfield(chain.params, 'ground')
    pi_mjp_ground =  (expm((chain.params.t_left(2)-chain.params.t_left(1))*squeeze(chain.params.ground.Q(:,:)))) ;
    end

    
    
end

pi11_mjp      = (squeeze(pi_mjp(1,1,:)));

pi22_mjp      = (squeeze(pi_mjp(2,2,:)));

pi33_mjp      = (squeeze(pi_mjp(3,3,:)));


axes1 = subplot(3,2,1);
p1    = plot(chain.i,pi11_mjp,...
    chain.i(idx),pi11_mjp(idx),'.-');
legend([ p1(2)],'{Escape Rate}','{True Escape Rate}');
if isfield(chain.params,'ground')
line(axes1.XLim,pi_mjp_ground(1,1).*[1 1],'linestyle','-','color','c','LineWidth',1);
end
ylabel('\pi_{1\rightarrow1} (1/s)');



axes11 = subplot(3,2,2);
histogram(pi11_mjp(idx),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
ylim(axes1.YLim);
if isfield(chain.params,'ground')
line(axes11.XLim,pi_mjp_ground(1,1).*[1 1],'linestyle','-','color','c','LineWidth',1);
end

axes2= subplot(3,2,3);
p1 = plot(chain.i,pi22_mjp,...
    chain.i(idx),pi22_mjp(idx),'.-');
 ylabel('k_{22} (1/s)');
if isfield(chain.params,'ground')
line(axes11.XLim,pi_mjp_ground(2,2).*[1 1],'linestyle','-','color','c','LineWidth',1);
end

axes21 = subplot(3,2,4);
histogram(pi22_mjp(idx),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
ylim(axes2.YLim);
if isfield(chain.params,'ground')
line(axes21.XLim,pi_mjp_ground(2,2).*[1 1],'linestyle','-','color','c','LineWidth',1);
end

axes3= subplot(3,2,5);
p1 = plot(chain.i,pi22_mjp,...
    chain.i(idx),pi22_mjp(idx),'.-');
 ylabel('k_{22} (1/s)');
if isfield(chain.params,'ground')
line(axes11.XLim,pi_mjp_ground(3,3).*[1 1],'linestyle','-','color','c','LineWidth',1);
end

axes31 = subplot(3,2,6);
histogram(pi22_mjp(idx),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
ylim(axes3.YLim);
if isfield(chain.params,'ground')
line(axes31.XLim,pi_mjp_ground(3,3).*[1 1],'linestyle','-','color','c','LineWidth',1);
end

