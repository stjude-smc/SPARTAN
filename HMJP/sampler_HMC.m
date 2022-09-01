function [mu_D,mu_A,mu_acc_rate_HMC,rec] = sampler_HMC(...
            ...
            ...
            t_s,t_f,...
            mu_D_old,mu_A_old,...
            mu_acc_rate_HMC,...
            iS,B,...
            params)
                                                
HMC_eps  = params.HMC_eps*rand(1);
HMC_L    = params.HMC_L;
demo     = false;

% MASSES
m_MSp_D = ones(1,B);
m_MSp_A = ones(1,B);

% MOMENTUM
sample_MSp_D = m_MSp_D.*randn(1,B);
sample_MSp_A = m_MSp_D.*randn(1,B);

% LEAP-FROG PREP
propos_mSq_D  = nan(HMC_L,B);
propos_mSq_A  = nan(HMC_L,B);

propos_MSp_D  = nan(HMC_L,B);
propos_MSp_A  = nan(HMC_L,B);

% FIRST STEP
i = 1;
mu_D_old          = mu_D_old';
propos_mSq_D(i,:) = mu_D_old;
mu_A_old          = mu_A_old';
propos_mSq_A(i,:) = mu_A_old;

params.dt = params.t_left(2)-params.t_left(1);

t_f = reshape( t_f,length(params.Int_D),1,length(t_s) );


   [U_mS_D,U_mS_A]= find_U_grad(mu_D_old,mu_A_old,t_s,t_f,iS,B,params);
propos_MSp_D(i,:) = sample_MSp_D - 0.5*HMC_eps*U_mS_D;
propos_MSp_A(i,:) = sample_MSp_A - 0.5*HMC_eps*U_mS_A;

if demo;Gim = HMC_visual([],B,HMC_L,propos_mSq_D,propos_mSq_A,propos_MSp_D,propos_MSp_A); end


% leap-frog forward
for i = 2:HMC_L-1
    propos_mSq_D(i,:) = propos_mSq_D(i-1,:) + HMC_eps*propos_MSp_D(i-1,:)./m_MSp_D;
    propos_mSq_A(i,:) = propos_mSq_A(i-1,:) + HMC_eps*propos_MSp_A(i-1,:)./m_MSp_A;
   
      [U_mS_D,U_mS_A]  = find_U_grad(propos_mSq_D(i,:),propos_mSq_A(i,:),t_s,t_f,iS,B,params);
    propos_MSp_D(i,:) = propos_MSp_D(i-1,:) - HMC_eps*U_mS_D;
    propos_MSp_A(i,:) = propos_MSp_A(i-1,:) - HMC_eps*U_mS_A;
    
    if demo;HMC_visual(Gim,B,[],propos_mSq_D,propos_mSq_A,propos_MSp_D,propos_MSp_A); end
end

% term step
i = HMC_L;
propos_mSq_D(i,:) = propos_mSq_D(i-1,:) + HMC_eps*propos_MSp_D(i-1,:)./m_MSp_D;
propos_mSq_A(i,:) = propos_mSq_A(i-1,:) + HMC_eps*propos_MSp_A(i-1,:)./m_MSp_A;

[U_mS_D,U_mS_A] = find_U_grad(propos_mSq_D(i,:),propos_mSq_A(i,:),t_s,t_f,iS,B,params);
propos_MSp_D(i,:) = propos_MSp_D(i-1,:) - 0.5*HMC_eps*U_mS_D;
propos_MSp_A(i,:) = propos_MSp_A(i-1,:) - 0.5*HMC_eps*U_mS_A;

if demo;HMC_visual(Gim,B,[],propos_mSq_D,propos_mSq_A,propos_MSp_D,propos_MSp_A); end


log_a = find_U(mu_D_old,mu_A_old,t_s,t_f,iS,params) ...
      - find_U(propos_mSq_D(end,:),propos_mSq_A(end,:),t_s,t_f,iS,params) ...
      + 0.5*( (sample_MSp_D.^2-propos_MSp_D(end,:).^2)./m_MSp_D...
             +(sample_MSp_A.^2-propos_MSp_A(end,:).^2)./m_MSp_A );
         
if log(rand) < log_a
    mu_D      = propos_mSq_D(end,:)';
    mu_A      = propos_mSq_A(end,:)';
    mu_acc_rate_HMC(1) = mu_acc_rate_HMC(1)+1;  
else
    mu_D = mu_D_old';
    mu_A = mu_A_old';
end
    mu_acc_rate_HMC(2) = mu_acc_rate_HMC(2)+1;
    
if ~iscolumn(mu_D)
mu_D = mu_D';
end

end

%% ploter
function Gim = HMC_visual(Gim,B,HMC_L,propos_mS_D_q,...
                                     propos_mS_A_q, ...
                                     propos_MS_D_p,...
                                     propos_MS_A_p)

if isempty(Gim)
    figure(88)
    
    col = [0 1 0;%green
           0 0 1;%blue
           1 0 0];%red
    
    subplot(1,3,[1 2])
    Gim.ax_Di{1} = plot(1:HMC_L,propos_mS_D_q,'o-'        ,'color',col(1,:));
    xlim([0 HMC_L+1])
    Gim.ax_Di{2} = line(1:HMC_L,propos_mS_A_q,'marker','o','color',col(3,:));
    xlim([0 HMC_L+1])
  
    subplot(1,3,3)
    Gim.ax_Dm{1} = plot(propos_MS_D_p,propos_mS_D_q,'o-'        ,'color',col(1,:));
    Gim.ax_Dm{2} = line(propos_MS_A_p,propos_mS_A_q,'marker','o'   ,'color',col(3,:));

end

for m=1:B
    Gim.ax_Di{1}(m).YData = propos_mS_D_q(:,m);
    Gim.ax_Dm{1}(m).XData = propos_MS_D_p(:,m);
    Gim.ax_Dm{1}(m).YData = propos_mS_D_q(:,m);
    Gim.ax_Di{2}(m).YData = propos_mS_A_q(:,m);
    Gim.ax_Dm{2}(m).XData = propos_MS_A_p(:,m);
    Gim.ax_Dm{2}(m).YData = propos_mS_A_q(:,m);
end

drawnow

end


%% HMC potential
function U = find_U(mS_D,mS_A,t_s,t_f,S,params)

if any(mS_D<0)||any(mS_A<0)
    U = inf;
else

       
       
    mS_D_temp  = reshape(mS_D(S),1,1,length(t_s)); % obs*1*length(t_t)
    mS_A_temp  = reshape(mS_A(S),1,1,length(t_s)); % obs*1*length(t_t)
    
    INT_MEAN_D  = (params.dt)*...
               (params.mu_back_D + sum(mS_D_temp.*t_f,3) );


    INT_MEAN_A  = (params.dt)*...
               (params.mu_back_A + sum(mS_A_temp.*t_f,3) );    
           
           
    U_D = -sum(  params.Int_D'.*log(INT_MEAN_D)-INT_MEAN_D) ...    
        + (1-params.phi_D)*sum( log(mS_D) ) + params.phi_D/params.psi_D*sum( mS_D );
    
    U_A = -sum(  params.Int_A'.*log(INT_MEAN_A)-INT_MEAN_A) ...    
            + (1-params.phi_A)*sum( log(mS_A) ) + params.phi_A/params.psi_A*sum( mS_A );
   
    U   = U_D+U_A;
end


end
    
%% HMC potential gradient
function [U_mS_D,U_mS_A] = find_U_grad(mS_D,mS_A,t_s,t_f,S,B,params)

      
  mS_D_temp  = reshape(mS_D(S),1,1,length(t_s)); % obs*1*length(t_t)
  mS_A_temp  = reshape(mS_A(S),1,1,length(t_s)); % obs*1*length(t_t)
         

U_mS_D = nan(1,B);

for m = 1:B
      ts_r    = reshape(S',1,1,length(t_s));

    U_mS_D(m) = -sum(( params.Int_D'.*(sum(( ts_r==m ).*t_f,3)))./(params.mu_back_D + (sum(mS_D_temp.*t_f,3)) )...
                - (params.dt)*(sum(( ts_r==m ).*t_f,3)))...
                + (1-params.phi_D)/mS_D(m) ...
                + params.phi_D/params.psi_D;
          
 
end

U_mS_A = nan(1,B);

for m = 1:B
    ts_r    = reshape(t_s,1,1,length(t_s));

    U_mS_A(m) = -sum(( params.Int_A'.*(sum(( ts_r==m ).*t_f,3)))./(params.mu_back_A + (sum(mS_A_temp.*t_f,3)) )...
                - (params.dt)*(sum(( ts_r==m ).*t_f,3)))...
                + (1-params.phi_A)/mS_A(m) ...
                + params.phi_A/params.psi_A;
   
     
end

end














% function [mu_D,mu_A,mu_acc_rate_HMC] = sampler_HMC(...
%             ...
%             ...
%             t_s,t_f,...
%             mu_D_old,mu_A_old,...
%             mu_acc_rate_HMC,...
%             params)
%                                                 
% HMC_eps  = params.HMC_eps*rand(1);
% HMC_L    = params.HMC_L;
% demo     = false;
% 
% % MASSES
% m_MSp_D = ones(1,params.M);
% m_MSp_A = ones(1,params.M);
% 
% % MOMENTUM
% sample_MSp_D = m_MSp_D.*randn(1,params.M);
% sample_MSp_A = m_MSp_D.*randn(1,params.M);
% 
% % LEAP-FROG PREP
% propos_mSq_D  = nan(HMC_L,params.M);
% propos_mSq_A  = nan(HMC_L,params.M);
% 
% propos_MSp_D  = nan(HMC_L,params.M);
% propos_MSp_A  = nan(HMC_L,params.M);
% 
% % FIRST STEP
% i = 1;
% mu_D_old          = mu_D_old';
% propos_mSq_D(i,:) = mu_D_old;
% mu_A_old          = mu_A_old';
% propos_mSq_A(i,:) = mu_A_old;
% 
% 
% 
%    [U_mS_D,U_mS_A]= find_U_grad(mu_D_old,mu_A_old,t_s,t_f,params);
% propos_MSp_D(i,:) = sample_MSp_D - 0.5*HMC_eps*U_mS_D;
% propos_MSp_A(i,:) = sample_MSp_A - 0.5*HMC_eps*U_mS_A;
% 
% if demo;Gim = HMC_visual([],params.M,HMC_L,propos_mSq_D,propos_mSq_A,propos_MSp_D,propos_MSp_A); end
% 
% 
% % leap-frog forward
% for i = 2:HMC_L-1
%     propos_mSq_D(i,:) = propos_mSq_D(i-1,:) + HMC_eps*propos_MSp_D(i-1,:)./m_MSp_D;
%     propos_mSq_A(i,:) = propos_mSq_A(i-1,:) + HMC_eps*propos_MSp_A(i-1,:)./m_MSp_A;
%    
%       [U_mS_D,U_mS_A]  = find_U_grad(propos_mSq_D(i,:),propos_mSq_A(i,:),t_s,t_f,params);
%     propos_MSp_D(i,:) = propos_MSp_D(i-1,:) - HMC_eps*U_mS_D;
%     propos_MSp_A(i,:) = propos_MSp_A(i-1,:) - HMC_eps*U_mS_A;
%     
%     if demo;HMC_visual(Gim,params.M,[],propos_mSq_D,propos_mSq_A,propos_MSp_D,propos_MSp_A); end
% end
% 
% % term step
% i = HMC_L;
% propos_mSq_D(i,:) = propos_mSq_D(i-1,:) + HMC_eps*propos_MSp_D(i-1,:)./m_MSp_D;
% propos_mSq_A(i,:) = propos_mSq_A(i-1,:) + HMC_eps*propos_MSp_A(i-1,:)./m_MSp_A;
% 
% [U_mS_D,U_mS_A] = find_U_grad(propos_mSq_D(i,:),propos_mSq_A(i,:),t_s,t_f,params);
% propos_MSp_D(i,:) = propos_MSp_D(i-1,:) - 0.5*HMC_eps*U_mS_D;
% propos_MSp_A(i,:) = propos_MSp_A(i-1,:) - 0.5*HMC_eps*U_mS_A;
% 
% if demo;HMC_visual(Gim,params.M,[],propos_mSq_D,propos_mSq_A,propos_MSp_D,propos_MSp_A); end
% 
% 
% log_a = find_U(mu_D_old,mu_A_old,t_s,t_f,params) ...
%       - find_U(propos_mSq_D(end,:),propos_mSq_A(end,:),t_s,t_f,params) ...
%       + 0.5*( (sample_MSp_D.^2-propos_MSp_D(end,:).^2)./m_MSp_D...
%              +(sample_MSp_A.^2-propos_MSp_A(end,:).^2)./m_MSp_A );
%          
% if log(rand) < log_a
%     mu_D      = propos_mSq_D(end,:)';
%     mu_A      = propos_mSq_A(end,:)';
%     mu_acc_rate_HMC(1) = mu_acc_rate_HMC(1)+1;  
% else
%     mu_D = mu_D_old';
%     mu_A = mu_A_old';
% end
%     mu_acc_rate_HMC(2) = mu_acc_rate_HMC(2)+1;
%     
% if ~iscolumn(mu_D)
% mu_D = mu_D';
% end
% 
% end
% 
% %% ploter
% function Gim = HMC_visual(Gim,M,HMC_L,propos_mS_D_q,...
%                                      propos_mS_A_q, ...
%                                      propos_MS_D_p,...
%                                      propos_MS_A_p)
% 
% if isempty(Gim)
%     figure(88)
%     
%     col = [0 1 0;%green
%            0 0 1;%blue
%            1 0 0];%red
%     
%     subplot(1,3,[1 2])
%     Gim.ax_Di{1} = plot(1:HMC_L,propos_mS_D_q,'o-'        ,'color',col(1,:));
%     xlim([0 HMC_L+1])
%     Gim.ax_Di{2} = line(1:HMC_L,propos_mS_A_q,'marker','o','color',col(3,:));
%     xlim([0 HMC_L+1])
%   
%     subplot(1,3,3)
%     Gim.ax_Dm{1} = plot(propos_MS_D_p,propos_mS_D_q,'o-'        ,'color',col(1,:));
%     Gim.ax_Dm{2} = line(propos_MS_A_p,propos_mS_A_q,'marker','o'   ,'color',col(3,:));
% 
% end
% 
% for m=1:M
%     Gim.ax_Di{1}(m).YData = propos_mS_D_q(:,m);
%     Gim.ax_Dm{1}(m).XData = propos_MS_D_p(:,m);
%     Gim.ax_Dm{1}(m).YData = propos_mS_D_q(:,m);
%     Gim.ax_Di{2}(m).YData = propos_mS_A_q(:,m);
%     Gim.ax_Dm{2}(m).XData = propos_MS_A_p(:,m);
%     Gim.ax_Dm{2}(m).YData = propos_mS_A_q(:,m);
% end
% 
% drawnow
% 
% end
% 
% 
% %% HMC potential
% function U = find_U(mS_D,mS_A,t_s,t_f,params)
% 
% if any(mS_D<0)||any(mS_A<0)
%     U = inf;
% else
% 
%        ts_r    = reshape(t_s,1,1,length(t_s));
%      
%     INT_MEAN_D = (params.t_right(1)-params.t_left(1))*...
%                (params.mu_back_D + mS_D(1)*(sum((ts_r==1).*t_f,3))+mS_D(2)*(1-(sum((ts_r==1).*t_f,3))));
% 
%     INT_MEAN_A = (params.t_right(1)-params.t_left(1))*...
%                (params.mu_back_A+ mS_A(1)*(sum((ts_r==1).*t_f,3))+mS_A(2)*(1-(sum((ts_r==1).*t_f,3))));
% 
%            
%     U_D = -sum( params.Int_D'.*log(INT_MEAN_D)-INT_MEAN_D) ...    
%         + (1-params.phi_D)*sum( log(mS_D) ) + params.phi_D/params.psi_D*sum( mS_D );
%     
%     U_A = -sum( params.Int_A'.*log(INT_MEAN_A)-INT_MEAN_A) ...    
%             + (1-params.phi_A)*sum( log(mS_A) ) + params.phi_A/params.psi_A*sum( mS_A );
%    
%     U   = U_D+U_A;
% end
% 
% 
% end
%     
% %% HMC potential gradient
% function [U_mS_D,U_mS_A] = find_U_grad(mS_D,mS_A,t_s,t_f,params)
% 
%       ts_r   = reshape(t_s,1,1,length(t_s));
%          
% INT_MEAN_A_d = params.Int_A'./(params.mu_back_A + mS_A(1)*(sum((ts_r==1).*t_f,3))+mS_A(2)*(1-(sum((ts_r==1).*t_f,3))) )...
%              -(params.t_right(1)-params.t_left(1));
% 
% U_mS_D = nan(1,params.M);
% 
% for m = 1:params.M
%       ts_r    = reshape(t_s,1,1,length(t_s));
% 
%     U_mS_D(m) = -sum((params.Int_D'.*(sum(( ts_r==m ).*t_f,3)))./(params.mu_back_D + mS_D(1)*(sum((ts_r==1).*t_f,3))+mS_D(2)*(1-(sum((ts_r==1).*t_f,3))) )...
%                  - (params.t_right(1)-params.t_left(1))*(sum(( ts_r==m ).*t_f,3)))...
%               + (1-params.phi_D)/mS_D(m) ...
%               + params.phi_D/params.psi_D;
% end
% 
% U_mS_A = nan(1,params.M);
% 
% for m = 1:params.M
%     ts_r    = reshape(t_s,1,1,length(t_s));
% 
%     U_mS_A(m) = -sum((params.Int_A'.*(sum(( ts_r==m ).*t_f,3)))./(params.mu_back_A + mS_A(1)*(sum((ts_r==1).*t_f,3))+mS_A(2)*(1-(sum((ts_r==1).*t_f,3))) )...
%                  - (params.t_right(1)-params.t_left(1))*(sum(( ts_r==m ).*t_f,3)))...
%               + (1-params.phi_A)/mS_A(m) ...
%               + params.phi_A/params.psi_A;
% end
% 
% end
% 
% 

