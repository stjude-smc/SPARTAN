function Gim = chainer_visualize(Gim,chain)



%% -------------------------------------------------------------------------
if isempty(Gim)

addpath('visualize')

figure(10)
set(gcf,'windowstyle','docked');


% LABELS AND AXES
% 
% t_true      = chain.params.ground.t_t;
% s_true      = chain.params.ground.t_s; 
% % 
% [tx_true,sx_true] = stairs(t_true,s_true);

t_t         =  (chain.sample.t_t);
t_s         =  (chain.sample.t_s);


    if chain.sample.t_t(end) < chain.params.T_f
        chain.sample.t_t = [chain.sample.t_t; chain.params.T_f];
        chain.sample.t_s = [chain.sample.t_s;chain.sample.t_s(end)];
    end
    
    

axes1 = subplot(1,1,1);
[tx_D,sx_D]     =  stairs( chain.sample.t_t,(chain.params.mu_back_D+chain.sample.mu_D(chain.sample.t_s))*...
                                        (chain.params.t_right(1)-chain.params.t_left(1)));
                                    
[tx_A,sx_A]     =  stairs( chain.sample.t_t,(chain.params.mu_back_A+chain.sample.mu_A(chain.sample.t_s))*...
                                        (chain.params.t_right(1)-chain.params.t_left(1)));
                                    
Gim.p1 = plot(tx_D,sx_D,'-o','color',[0.8 0.8 0.8],'LineWidth',2);
hold on
Gim.p2 = plot(tx_A,sx_A,'-o','color','c','LineWidth',2);

line(chain.params.T_i.*[1 1],[0,max(chain.sample.mu_D(chain.sample.t_s))+0.1],'linestyle','--','color','k','linewidth',2);
line(chain.params.T_f.*[1 1],[0,max(chain.sample.mu_D(chain.sample.t_s))+0.1],'linestyle','--','color','k','linewidth',2);

% l1 = line(tx_true,(chain.params.mu_back_D+chain.params.ground.mu_D(sx_true))*...
%              (chain.params.t_right(1)-chain.params.t_left(1)),'linestyle','-','color','g')
%  
% l2 = line(tx_true,(chain.params.mu_back_A+chain.params.ground.mu_A(sx_true))*...
%              (chain.params.t_right(1)-chain.params.t_left(1)),'linestyle','-','color','r')
        
% line([chain.params.T_i,chain.params.T_f],chain.params.ground.mu_D([1;2])'.*[1 ;1],'linestyle','-','color','c','LineWidth',.1);

    for j=1:length(chain.params.t_left)
    
        line([chain.params.t_left(j) chain.params.t_right(j)],...
            [chain.params.Int_D(j) chain.params.Int_D(j)],'color',[0 .7 0],'LineWidth',5);
    
        line([chain.params.t_left(j) chain.params.t_right(j)],...
            [chain.params.Int_A(j) chain.params.Int_A(j)],'color',[.7 0 0],'LineWidth',5);
    end
        hold off
    leg = legend([Gim.p1(1) Gim.p2(1) ],[{'Donor (MCMC)','Acceptor (MCMC)'}, {'Donor (True)','Acceptor (True)'}],'location','NorthEast','orientation','horizontal')
    legend boxoff
    leg.NumColumns = 2;
    leg.Position = [0.064540471732356,0.946707824908572,0.455160736079547,0.035185184302153];
xlim_0   =  chain.params.T_i - 0.05*(chain.params.T_f-chain.params.T_i);
xlim_end =  chain.params.T_f + 0.05*(chain.params.T_f-chain.params.T_i) ;


xlim([xlim_0 xlim_end]);
ylim([0,max([(chain.params.mu_back_D(end)+chain.sample.mu_D([1;2],end))*(-chain.params.t_left(1)+chain.params.t_right(1));...
    (chain.sample.mu_A([1;2],(end))+chain.params.mu_back_A((end)))*(-chain.params.t_left(1)+chain.params.t_right(1))])+50]);

xlabel('Time (s)');
ylabel('Intensity\newline(photons)');

set(gca);

else

t_t         =  (chain.sample.t_t);
t_s         =  (chain.sample.t_s);


 if chain.sample.t_t(end) < chain.params.T_f
          chain.sample.t_t = [chain.sample.t_t; chain.params.T_f];
          chain.sample.t_s = [chain.sample.t_s;chain.sample.t_s(end)];
 end

[tx_D,sx_D]     =  stairs( chain.sample.t_t,(chain.params.mu_back_D+chain.sample.mu_D(chain.sample.t_s))*...
                                        (chain.params.t_right(1)-chain.params.t_left(1)) );    
[tx_A,sx_A]     =  stairs( chain.sample.t_t,(chain.params.mu_back_A+chain.sample.mu_A(chain.sample.t_s))*...
                                        (chain.params.t_right(1)-chain.params.t_left(1)) );    
    
set(Gim.p1,'XData',tx_D,'YData',sx_D);
hold on
set(Gim.p2,'XData',tx_A,'YData',sx_A);

ylimits = get(gca,'yLim');


yvals   = linspace(ylimits(1),ylimits(end),100);
mu_D_prior = 1/((chain.params.psi_D/chain.params.phi_D)^chain.params.phi_D)*...
              (gamma(chain.params.phi_D))*((yvals).^(chain.params.phi_D-1)).*...
              exp(-(yvals*chain.params.phi_D/chain.params.psi_D));


        
hold off
end



drawnow