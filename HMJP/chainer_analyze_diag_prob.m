function [pi11_mjp ,pi22_mjp] = chainer_analyze_diag_prob(chain,fr,dr,m_num,demo_flag)
%  load('')

%% collect samples
r_max = chain.length;
r_min = fix( fr*r_max );
idx = r_min:dr:r_max;


pi_mjp        = nan(chain.params.M,chain.params.M,(chain.length));

% pi_ground     = (expm((chain.params.t_left(2)-chain.params.t_left(1))*chain.params.ground.Q));
for j=1:(chain.length)   
    
pi_mjp(:,:,j) = (expm((chain.params.t_left(2)-chain.params.t_left(1))*squeeze(chain.Q(:,:,j)))) ;

end

pi11_mjp      = (squeeze(pi_mjp(1,1,:)));

pi22_mjp      = (squeeze(pi_mjp(2,2,:)));

if demo_flag
    
  
figure


axes11 = subplot(2,2,1);
p1    = plot(chain.i,pi11_mjp,...
    chain.i(idx),pi11_mjp(idx),'.-');
% l1    = line(axes1.XLim,pi_ground(1,1).*[1 1],'linestyle','-','color','c','LineWidth',1);
legend([ p1(2)],'{Diagonal transition probability}');

ylabel('\pi_{11} ');



axes12 = subplot(2,2,2);
histogram(pi11_mjp(idx),m_num,'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
ylim(axes11.YLim);
% line(axes11.XLim,pi_ground(1,1).*[1 1],'linestyle','-','color','c','LineWidth',1);


axes21= subplot(2,2,3);
p1 = plot(chain.i,pi22_mjp,...
    chain.i(idx),pi22_mjp(idx),'.-');
 ylabel('\pi_{22} ');


axes22 = subplot(2,2,4);
histogram(pi22_mjp(idx),m_num,'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
ylim(axes22.YLim);

end

% axes3 = subplot(3,2,5);
% p1 = plot(chain.i,pi33_mjp,...
%     chain.i(idx),pi33_mjp(idx),'.-');
%     line(axes3.XLim,pi_ground(3,3).*[1 1],'linestyle','-','color','c','LineWidth',1);
% ylabel('k_{33} (1/s)');
% 
% xlabel('MCMC Chain')
% 
% axes31 = subplot(3,2,6);
% histogram(pi33_mjp(idx),'normalization','pdf','Orientation','horizontal','Facecolor',get(p1(2),'color'));
% ylim(axes3.YLim);
% line(axes31.XLim,pi_ground(3,3).*[1 1],'linestyle','-','color','c','LineWidth',1);
% xlabel('PDF');

















%% Ground Truth
% 
% Q_g             = - [-3     0.9  2.1 ;...
%                    3.8  -6.1  2.3;...
%                    3.94  2.06 -6 ];

%% Visualizing escape rates and Mean Dwell Times
% 
% s     = size(chain.Q)
% 
% for j=2:s(3)
%     lambda1(j-1) = -chain.Q(1,1,j);
% end
% for j=1001:s(3)
%     lambda2(j-1000) = -chain.Q(2,2,j);
% end
% for j=1001:s(3)
%     lambda3(j-1000) = -chain.Q(3,3,j);
% end
% 
% figure(55)
% title('Escape Rates and Mean Dwell Times')
% 
% maxrate = max([lambda1 lambda2 lambda3]);
% minrate = min([lambda1 lambda2 lambda3]);
% 
% axes11 = subplot(3,2,1)
% hold on
% title('Escape Rates')
% histogram(lambda1,linspace(0,maxrate,100),'normalization','pdf');
% yLim = axes11.YLim;
% line(median(lambda1).*[1 1],yLim,'linestyle','--','color','r','linewidth',2)
% line(Q_g(1,1).*[1 1],yLim,'linestyle',':','color','b','linewidth',2)
% ytemp = linspace(0.1,5,100);
% escrate_prior1 = plot(ytemp,gampdf(ytemp,4/2,2/(4*1)),'linestyle','--','color','m','linewidth',2)
% xlabel('\lambda_1')
% xlim([0,10])
% ylabel('PDF')
% hold off
% legend('learned','median','true','prior','location','NE')
% 
% axes12 = subplot(3,2,2)
% hold on
% title('Mean Dwell Times')
% histogram(1./lambda1,linspace(0,1/minrate,100),'normalization','pdf');
% line(median(1./lambda1).*[1 1],axes12.YLim,'linestyle','--','color','r','linewidth',2)
% line((1/Q_g(1,1))*[1 1],axes12.YLim,'linestyle',':','color','b','linewidth',2)
% 
% xlabel('1/\lambda_1')
% set(axes12,'XScale','lin')
% ylabel('PDF')
% hold off
% legend('learned','median','true','location','NE')
% 
% %%%%%%%
% axes21 = subplot(3,2,3)
% histogram(lambda2,linspace(0,maxrate,100),'normalization','pdf');
% line(median(lambda2).*[1 1],axes21.YLim,'linestyle','--','color','r','linewidth',2)
% line(Q_g(2,2).*[1 1],axes21.YLim,'linestyle',':','color','b','linewidth',2)
% ytemp = linspace(0.1,5,100);
% hold on
% escrate_prior2 = plot(ytemp,gampdf(ytemp,4/2,2/(4*1)),'linestyle','--','color','m','linewidth',2)
% xlabel('\lambda_2')
% ylabel('PDF')
% legend('learned','median','true','prior','location','NE')
% 
% axes22 = subplot(3,2,4)
% histogram(1./lambda2,linspace(0,1/minrate,100),'normalization','pdf');
% line(median(1./lambda2).*[1 1],axes22.YLim,'linestyle','--','color','r','linewidth',2)
% line((1./(Q_g(2,2))).*[1 1],axes22.YLim,'linestyle',':','color','b','linewidth',2)
% set(axes22,'XScale','lin')
% 
% xlabel('1/\lambda_2')
% ylabel('PDF')
% hold off
% legend('learned','median','true','location','NE')
% 
% axes31 = subplot(3,2,5)
% histogram(lambda3,linspace(0,maxrate,100),'normalization','pdf');
% line(median(lambda3).*[1 1],axes31.YLim,'linestyle','--','color','r','linewidth',2)
% line(Q_g(3,3).*[1 1],axes31.YLim,'linestyle',':','color','b','linewidth',2)
% ytemp = linspace(0.1,5,100);
% hold on
% escrate_prior = plot(ytemp,gampdf(ytemp,4/2,2/(4*1)),'linestyle','--','color','m','linewidth',2)
% xlabel('\lambda_3')
% ylabel('PDF')
% hold off
% % legend('learned','median','true','prior','location','NE')
% 
% axes32 = subplot(3,2,6)
% histogram(1./lambda3,linspace(0,1/minrate,100),'normalization','pdf');
% line(median(1./lambda3).*[1 1],axes32.YLim,'linestyle','--','color','r','linewidth',2)
% line(median(1./Q_g(3,3)).*[1 1],axes32.YLim,'linestyle',':','color','b','linewidth',2)
% set(axes32,'XScale','lin')
% 
% xlabel('1/\lambda_3')
% ylabel('PDF')
% legend('learned','median','true','location','NE')
