function [mu_D,mu_A,mu_acc_rate_MH,mu_acc_rate_HMC,acc_rate_flip] = sampler_update_emission(...
    ...
    ...
    t_s,t_f,...
    mu_D,mu_A,...
    mu_acc_rate_MH,...
    mu_acc_rate_HMC,...
    acc_rate_flip,...
    params)
% this function updates the molecule state levels


[S,~,iS] = unique(t_s);    
 B = length(S);

% Propose mean_st from a gamma distribution
for r = 1:5
        [mu_D(S),mu_A(S),mu_acc_rate_HMC] = sampler_HMC(...
                                                  ...
                                                  ...
                                                  t_s,t_f,...
                                                  mu_D(S),mu_A(S),...
                                                  mu_acc_rate_HMC,...
                                                  iS,B,...
                                                  params);
        [mu_D(S),mu_A(S),mu_acc_rate_MH] = sampler_MH(...
                                                ...
                                                ...
                                                t_s,t_f,...
                                                mu_D(S),mu_A(S),...
                                                mu_acc_rate_MH,...
                                                iS,...
                                                params);
        [mu_D(S),mu_A(S),acc_rate_flip] = sampler_flip(...
                                                ...
                                                ...
                                                t_s,t_f,...
                                                mu_D(S),mu_A(S),...
                                                acc_rate_flip,...
                                                iS,...
                                                params);
                                                                                            
                                                
                                                
end

% update inactive emission levels

S  = setdiff(1:params.M,S);
B = length(S);

mu_D(S)     = (params.psi_D/params.phi_D)*randg(params.phi_D,[B,1]);
mu_A(S)     = (params.psi_A/params.phi_A)*randg(params.phi_A,[B,1]);

end




function [mu_D,mu_A,mu_acc_rate] = sampler_MH(...
    ...
    ...
    t_s,t_f,...
    mu_D_old,mu_A_old,mu_acc_rate,...
    S,...
    params)

params.dt = params.t_left(2)-params.t_left(1);
% this function updates the molecule state levels
% Propose mean_st from a gamma distribution
for rr= 1:5
    switch randi(1)
        case 1
            propos_mu_D        = (mu_D_old/params.alpha_prop_D).*randg(params.alpha_prop_D,[length(mu_D_old) 1])    ; 
            propos_mu_A        = (mu_A_old/params.alpha_prop_A).*randg(params.alpha_prop_A,[length(mu_D_old) 1])    ; 
      
    end

%             propos_mean_st(1,1)   =  mean_st_old(1); 
%             propos_mean_st(2,1)   =  mean_st_old(2);
%     
 % Calculate the likelihood
 mu_D_s_temp_new            = reshape(propos_mu_D(S)/params.mu_back_D,1,1,length(S)); % obs*1*length(t_t)
 mu_D_s_temp_old            = reshape(mu_D_old(S)/params.mu_back_D,1,1,length(S)); % obs*1*length(t_t)
 
 mu_A_s_temp_new            = reshape(propos_mu_A(S)/params.mu_back_A,1,1,length(S)); % obs*1*length(t_t)
 mu_A_s_temp_old            = reshape(mu_A_old(S)/params.mu_back_A,1,1,length(S)); % obs*1*length(t_t)
    constant_like_D = 0;
    constant_like_A = 0;
   
                 neww_D      = (1+...
                              sum(mu_D_s_temp_new.*t_f,3)) ;
                          
                 oldd_D      = (1+...
                              sum(mu_D_s_temp_old.*t_f,3)) ;
                          
                 neww_A      = (1+...
                              sum(mu_A_s_temp_new.*t_f,3)) ;
                          
                 oldd_A      = (1+...
                              sum(mu_A_s_temp_old.*t_f,3)) ;         
                          
                 constant_like_D = sum( params.Int_D'.*(log(neww_D)-log(oldd_D))...
                                        +(params.mu_back_D*(params.dt)*(oldd_D-neww_D)));
                                
                 constant_like_A = sum( params.Int_A'.*(log(neww_A)-log(oldd_A))...
                                        +(params.mu_back_A*(params.dt)*(oldd_A-neww_A)));                
         if all(propos_mu_D>0) && all(propos_mu_A>0) 

         logr = constant_like_D + constant_like_A ...
                 ...
                +((2*params.alpha_prop_D-params.phi_D)*sum(log(mu_D_old./propos_mu_D)))...
                +(sum(mu_D_old-propos_mu_D)*(params.phi_D/params.psi_D))...
                +(params.alpha_prop_D*sum((propos_mu_D./mu_D_old)-(mu_D_old./propos_mu_D)))...
                ...
                +((2*params.alpha_prop_A-params.phi_A)*sum(log(mu_A_old./propos_mu_A)))...
                +(sum(mu_A_old-propos_mu_A)*(params.phi_A/params.psi_A))...
                +(params.alpha_prop_A*sum((propos_mu_A./mu_A_old)-(mu_A_old./propos_mu_A)));
        else
            logr = -inf;
        end
% Accept or reject the proposals
     if  logr>log(rand)
         mu_D             = propos_mu_D                ;
         mu_A             = propos_mu_A                ;
         mu_acc_rate(1)   = mu_acc_rate(1)+1           ;
     else
         mu_D             = mu_D_old                   ;
         mu_A             = mu_A_old                   ;
     end    
     mu_acc_rate(2)       = mu_acc_rate(2)+1           ; 
end  
end    


function [mu_D,mu_A,acc_rate_flip] = sampler_flip(...
    ...
    ...
    t_s,t_f,...
    mu_D_old,mu_A_old,acc_rate_flip,...
    S,...
    params)
% this function updates the molecule state levels

params.dt = params.t_left(2)-params.t_left(1);
B = length(unique(S));
if B>1
m1 = randi(B);
m2 = randi(B-1);
m2 = m2 + (m2>=m1);

propos_mu_D     = mu_D_old;
propos_mu_D(m2) = mu_D_old(m1);
propos_mu_D(m1) = mu_D_old(m2);


propos_mu_A = mu_A_old;
propos_mu_A(m2) = mu_A_old(m1);
propos_mu_A(m1) = mu_A_old(m2);

%     
 % Calculate the likelihood
 mu_D_s_temp_new            = reshape(propos_mu_D(S)/params.mu_back_D,1,1,length(S)); % obs*1*length(t_t)
 mu_D_s_temp_old            = reshape(mu_D_old(S)/params.mu_back_D,1,1,length(S)); % obs*1*length(t_t)
 
 mu_A_s_temp_new            = reshape(propos_mu_A(S)/params.mu_back_A,1,1,length(S)); % obs*1*length(t_t)
 mu_A_s_temp_old            = reshape(mu_A_old(S)/params.mu_back_A,1,1,length(S)); % obs*1*length(t_t)
    constant_like_D = 0;
    constant_like_A = 0;
   
                 neww_D      = (1+...
                              sum(mu_D_s_temp_new.*t_f,3)) ;
                          
                 oldd_D      = (1+...
                              sum(mu_D_s_temp_old.*t_f,3)) ;
                          
                 neww_A      = (1+...
                              sum(mu_A_s_temp_new.*t_f,3)) ;
                          
                 oldd_A      = (1+...
                              sum(mu_A_s_temp_old.*t_f,3)) ;         
                          
                 constant_like_D = sum( params.Int_D'.*(log(neww_D)-log(oldd_D))...
                                        +(params.mu_back_D*(params.dt)*(oldd_D-neww_D)));
                                
                 constant_like_A = sum( params.Int_A'.*(log(neww_A)-log(oldd_A))...
                                        +(params.mu_back_A*(params.dt)*(oldd_A-neww_A)));                
         if all(propos_mu_D>0) && all(propos_mu_A>0) 

         logr = constant_like_D + constant_like_A ...
                 ...
                +((2*params.alpha_prop_D-params.phi_D)*sum(log(mu_D_old./propos_mu_D)))...
                +(sum(mu_D_old-propos_mu_D)*(params.phi_D/params.psi_D))...
                +(params.alpha_prop_D*sum((propos_mu_D./mu_D_old)-(mu_D_old./propos_mu_D)))...
                ...
                +((2*params.alpha_prop_A-params.phi_A)*sum(log(mu_A_old./propos_mu_A)))...
                +(sum(mu_A_old-propos_mu_A)*(params.phi_A/params.psi_A))...
                +(params.alpha_prop_A*sum((propos_mu_A./mu_A_old)-(mu_A_old./propos_mu_A)));
        else
            logr = -inf;
        end
% Accept or reject the proposals
     if  logr>log(rand)
         mu_D             = propos_mu_D                ;
         mu_A             = propos_mu_A                ;
         acc_rate_flip(1)   = acc_rate_flip(1)+1           ;
     else
         mu_D             = mu_D_old                   ;
         mu_A             = mu_A_old                   ;
     end    
     acc_rate_flip(2)       = acc_rate_flip(2)+1           ; 
else
    display('only one active levels');
         mu_D             = mu_D_old                   ;
         mu_A             = mu_A_old                   ;
   acc_rate_flip(2)       = acc_rate_flip(2)+1           ; 

end    

end



















% function [mu_D,mu_A,mu_acc_rate_MH,mu_acc_rate_HMC] = sampler_update_emission(...
%     ...
%     ...
%     t_s,t_f,...
%     mu_D,mu_A,...
%     mu_acc_rate_MH,...
%     mu_acc_rate_HMC,...
%     params)
% % this function updates the molecule state levels
% 
% if isempty(t_s)
%     propos_mu_D = nan(params.M,1);
%     propos_mu_A = nan(params.M,1);
% 
%     propos_mu_D      = (params.psi_D/params.phi_D)*randg(params.phi_D,[params.M,1]);
%     propos_mu_A      = (params.psi_A/params.phi_A)*randg(params.phi_A,[params.M,1]);
% else
% 
% % Propose mean_st from a gamma distribution
%     for r = 1:5
%         [mu_D,mu_A,mu_acc_rate_HMC] = sampler_HMC(...
%                                                   ...
%                                                   ...
%                                                   t_s,t_f,...
%                                                   mu_D,mu_A,...
%                                                   mu_acc_rate_HMC,...
%                                                   params);
%         [mu_D,mu_A,mu_acc_rate_MH] = sampler_MH(...
%                                                 ...
%                                                 ...
%                                                 t_s,t_f,...
%                                                 mu_D,mu_A,...
%                                                 mu_acc_rate_MH,...
%                                                 params);
%                                                 
%                                                 
%     end
% end
% 
% end
% 
% 
% 
% function [mu_D,mu_A,mu_acc_rate] = sampler_MH(...
%     ...
%     ...
%     t_s,t_f,...
%     mu_D_old,mu_A_old,mu_acc_rate,...
%     params)
% % this function updates the molecule state levels
% 
% % Propose mean_st from a gamma distribution
%     switch randi(7)
%         case 1
%             propos_mu_D        = (mu_D_old/params.alpha_prop_D).*randg(params.alpha_prop_D,[params.M 1])    ; 
%             propos_mu_A        = (mu_A_old/params.alpha_prop_A).*randg(params.alpha_prop_A,[params.M 1])    ; 
%         case 2
%             propos_mu_D(1,1)   = mu_D_old(1,1);
%             propos_mu_D(2,1)   = (mu_D_old(2,1)/params.alpha_prop_D).*randg(params.alpha_prop_D,1)          ; 
%             
%             propos_mu_A(1,1)   =  mu_A_old(1,1);
%             propos_mu_A(2,1)   = (mu_A_old(2,1)/params.alpha_prop_A).*randg(params.alpha_prop_A,1)          ; 
%         case 3
%             propos_mu_D(1,1)   = (mu_D_old(1,1)/params.alpha_prop_D).*randg(params.alpha_prop_D,1)          ; 
%             propos_mu_D(2,1)   =  mu_D_old(2,1);
%             
%             propos_mu_A(1,1)   = (mu_A_old(1,1)/params.alpha_prop_A).*randg(params.alpha_prop_A,1)          ; 
%             propos_mu_A(2,1)   =  mu_A_old(2,1);
%         case 4
%             propos_mu_D(1,1)   = (mu_D_old(1,1)/params.alpha_prop_D).*randg(params.alpha_prop_D,1)          ; 
%             propos_mu_D(2,1)   =  mu_D_old(2,1);
%             
%             propos_mu_A        = (mu_A_old/params.alpha_prop_A).*randg(params.alpha_prop_A,[params.M 1])    ;             
%         case 5
%             propos_mu_D        = (mu_D_old/params.alpha_prop_D).*randg(params.alpha_prop_D,[params.M 1])    ; 
%             
%             propos_mu_A(1,1)   = (mu_A_old(1,1)/params.alpha_prop_A).*randg(params.alpha_prop_A,1)          ; 
%             propos_mu_A(2,1)   =  mu_A_old(2,1);
%        case 6
%             propos_mu_D(1,1)   =  mu_D_old(1,1);
%             propos_mu_D(2,1)   = (mu_D_old(2,1)/params.alpha_prop_D).*randg(params.alpha_prop_D,1)          ; 
%         
%             propos_mu_A        = (mu_A_old/params.alpha_prop_A).*randg(params.alpha_prop_A,[params.M 1])    ; 
%             
%        case 7
%             propos_mu_D        = (mu_D_old/params.alpha_prop_D).*randg(params.alpha_prop_D,[params.M 1])    ; 
%     
%             propos_mu_A(1,1)   =  mu_A_old(1,1); 
%             propos_mu_A(2,1)   = (mu_A_old(2,1)/params.alpha_prop_A).*randg(params.alpha_prop_A,1)          ; 
% 
%     end
% %             propos_mean_st(1,1)   =  mean_st_old(1); 
% %             propos_mean_st(2,1)   =  mean_st_old(2);
% %     
%  % Calculate the likelihood
%  mu_D_s_temp_new            = reshape(propos_mu_D(t_s)/params.mu_back_D,1,1,length(t_s)); % obs*1*length(t_t)
%  mu_D_s_temp_old            = reshape(mu_D_old(t_s)/params.mu_back_D,1,1,length(t_s)); % obs*1*length(t_t)
%  
%  mu_A_s_temp_new            = reshape(propos_mu_A(t_s)/params.mu_back_A,1,1,length(t_s)); % obs*1*length(t_t)
%  mu_A_s_temp_old            = reshape(mu_A_old(t_s)/params.mu_back_A,1,1,length(t_s)); % obs*1*length(t_t)
%     constant_like_D = 0;
%     constant_like_A = 0;
%    
%                  neww_D      = (1+...
%                               sum(mu_D_s_temp_new.*t_f,3)) ;
%                           
%                  oldd_D      = (1+...
%                               sum(mu_D_s_temp_old.*t_f,3)) ;
%                           
%                  neww_A      = (1+...
%                               sum(mu_A_s_temp_new.*t_f,3)) ;
%                           
%                  oldd_A      = (1+...
%                               sum(mu_A_s_temp_old.*t_f,3)) ;         
%                           
%                  constant_like_D = sum(params.Int_D'.*(log(neww_D)-log(oldd_D))...
%                                         +(params.mu_back_D*(params.t_right(1)-params.t_left(1))*(oldd_D-neww_D)));
%                                 
%                  constant_like_A = sum(params.Int_A'.*(log(neww_A)-log(oldd_A))...
%                                         +(params.mu_back_A*(params.t_right(1)-params.t_left(1))*(oldd_A-neww_A)));                
%          if all(propos_mu_D>0) && all(propos_mu_A>0) 
% 
%          logr = constant_like_D + constant_like_A ...
%                  ...
%                 +((2*params.alpha_prop_D-params.phi_D)*sum(log(mu_D_old./propos_mu_D)))...
%                 +(sum(mu_D_old-propos_mu_D)*(params.phi_D/params.psi_D))...
%                 +(params.alpha_prop_D*sum((propos_mu_D./mu_D_old)-(mu_D_old./propos_mu_D)))...
%                 ...
%                 +((2*params.alpha_prop_A-params.phi_A)*sum(log(mu_A_old./propos_mu_A)))...
%                 +(sum(mu_A_old-propos_mu_A)*(params.phi_A/params.psi_A))...
%                 +(params.alpha_prop_A*sum((propos_mu_A./mu_A_old)-(mu_A_old./propos_mu_A)));
%         else
%             logr = -inf;
%         end
% % Accept or reject the proposals
%      if  logr>log(rand)
%          mu_D             = propos_mu_D                ;
%          mu_A             = propos_mu_A                ;
%          mu_acc_rate(1)   = mu_acc_rate(1)+1           ;
%      else
%          mu_D             = mu_D_old                   ;
%          mu_A             = mu_A_old                   ;
%      end    
%      mu_acc_rate(2)       = mu_acc_rate(2)+1           ; 
% end  
%     
%     