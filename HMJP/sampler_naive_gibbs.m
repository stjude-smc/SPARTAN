function [t_t,t_s] = sampler_naive_gibbs(t_t,t_s,mu_D,mu_A,mu_back_D,mu_back_A,IT,B,t_f,params)


%% -------------------------ALL PANELS----------------------------------
% k is for panels

%% WHICH OBSERVATION WINDOW DOES CONTAIN THE k-th HOLDING STATE (t_k-t_k+1) ?
dt = params.t_right(1)-params.t_left(1);

for k = 1:length(t_t)%randperm(length(t_t)) 

    if abs(k-1)<eps          
        P_in       = IT';     %column matrix     
    else              
        P_in       = B(t_s(k-1),:)';
    end
    
    if abs(k-length(t_t))>eps                        
        P_out      = B(:,t_s(k+1));         
    else              
        P_out      = ones(params.M,1);
    end
   
    mu_D_s_temp        = repmat(reshape(mu_D(t_s),1,1,length(t_s)),1,params.M,1); % obs*1*length(t_t)
    mu_D_s_temp(1,:,k) = mu_D((1:1:params.M));                                    % only the k th step,no Markovianity
    log_ld_s_D           = sum(-dt*(mu_back_D+sum(mu_D_s_temp.*t_f,3))+...
                               params.Int_D'.*log(dt*(mu_back_D+sum(mu_D_s_temp.*t_f,3)))-...
                               gammaln(params.Int_D'+1));   % 1 by M
                           
    mu_A_s_temp        = repmat(reshape(mu_A(t_s),1,1,length(t_s)),1,params.M,1); % obs*1*length(t_t)
    mu_A_s_temp(1,:,k) = mu_A((1:1:params.M));                                    % only the k th step,no Markovianity
    log_ld_s_A         = sum(-dt*(mu_back_A+sum(mu_A_s_temp.*t_f,3))+...
                               params.Int_A'.*log(dt*(mu_back_A+sum(mu_A_s_temp.*t_f,3)))-...
                               gammaln(params.Int_A'+1));   % 1 by M                       
                           

    t_s(k)              = gumble_sample( log_ld_s_D'+log_ld_s_A'+log(P_in)+log(P_out));

 
end





uu= 9;
 
   
