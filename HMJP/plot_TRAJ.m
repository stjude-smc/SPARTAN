Figures = figure(55)
set(gcf,'windowstyle','docked');

% cd .. 
% cd data
% chain = load('tag_mean_newpriordist_0pt05_5_96001.mat')
% cd ..
% cd visualize
s     = size(chain.IT);
chain.length =  s(1);
dr = 0.1; 
idx = (floor(dr*chain.length):1:chain.length);



% LABELS AND AXES
t_true      = chain.params.ground.t_t;
s_true      = chain.params.ground.t_s; 
% 
[tx_true,sx_true] = stairs(t_true,s_true);

axes1 = subplot(1,2,[1 2]);
j = find(max(chain.MAP(idx,1))==chain.MAP(idx,1));
% for j = 1:length(idx)
[tx_D,sx_D]     =  stairs( chain.t_t{idx(j)},(chain.params.mu_back_D+chain.mu_D((chain.t_s{idx(j)}),j))*...
                                        (chain.params.t_right(1)-chain.params.t_left(1)));
                                    
[tx_A,sx_A]     =  stairs( chain.t_t{idx(j)},(chain.params.mu_back_A+chain.mu_A((chain.t_s{idx(j)}),j))*...
                                        (chain.params.t_right(1)-chain.params.t_left(1)));                           
Gim.p1 = plot(tx_D,sx_D,'.-','color','g','LineWidth',2);
hold on
Gim.p2 = plot(tx_A,sx_A,'.-','color','r','LineWidth',2);
l1 = line(tx_true,(chain.params.mu_back_D+chain.params.ground.mu_D(sx_true))*...
             (chain.params.t_right(1)-chain.params.t_left(1)),'linestyle','-','color','c')
 
l2 = line(tx_true,(chain.params.mu_back_A+chain.params.ground.mu_A(sx_true))*...
             (chain.params.t_right(1)-chain.params.t_left(1)),'linestyle','-','color','c')


% end
        
    for j=1:length(chain.params.t_left)
    
        line([chain.params.t_left(j) chain.params.t_right(j)],...
            [chain.params.Int_D(j) chain.params.Int_D(j)],'color',[0 0.7 0],'LineWidth',5);
    
        line([chain.params.t_left(j) chain.params.t_right(j)],...
            [chain.params.Int_A(j) chain.params.Int_A(j)],'color',[0.7 0 0],'LineWidth',5);
    end
        hold off
        
xlim_0   =  chain.params.T_i - 0.05*(chain.params.T_f-chain.params.T_i);
xlim_end =  chain.params.T_f + 0.05*(chain.params.T_f-chain.params.T_i) ;


xlim([xlim_0 xlim_end]);
ylim([0,max(chain.params.Int_D)]);

xlabel('Time (s)');
ylabel('Intensity\newline(photons)');
leg = legend([Gim.p1(1) Gim.p2(1) l1(1) l2(1)],[{'Donor MAP traj.','Acceptor MAP traj.'},{'Donor true traj.','Acceptor true traj.'}],'location','NorthEast','orientation','horizontal')
    legend boxoff
    leg.NumColumns = 2;
    leg.Position = [0.207074116342842,0.858695685303109,0.562616811344557,0.035185184302153];



set(gca);
