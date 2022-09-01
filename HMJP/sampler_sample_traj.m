function [t_t,t_s] =sampler_sample_traj(t_t,t_s,...
                    mu_D,mu_A,mu_back_D,mu_back_A,...
                    IT,B,t_f,...
                   params)




  
 [t_t,t_s] = sampler_naive_gibbs(t_t,t_s,mu_D,mu_A,mu_back_D,mu_back_A,IT,B,t_f,params);    
% else    
%   keyboard
%  [t_t,t_s] = chainer_FFBS(t_t,t_s,mu_D,mu_A,t_f,params,P_gr,group_st,n_st,l_obs);

% 
% for m = 1:length(idx)
%     
%     if length(idx{m})>=1
% 
%        flag
%         
%        [t_t,t_s] = chainer_naive_gibbs(t_t,t_s,IT,B,t_d,idx,params);
%        
%     else
%         
%        [t_t,t_s] = chainer_FFBS(t_t,t_s,IT,B,t_d,idx,params);
%        
%     end
%     
% end
% 


