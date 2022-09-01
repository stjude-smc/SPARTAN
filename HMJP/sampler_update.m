function sample = sampler_update( sample, params )

%% update counter
sample.i = sample.i + 1;

%% Uniformize
[B_prob,omega]                          = sampler_uniformize(sample.Q,params);

%% Add Virtual Jumps
[sample.t_t, sample.t_s]                = sampler_add_jumps(sample.t_t,sample.t_s,sample.Q,omega,params);

%% Holding Times and Indexes AFTER the Virtual Jumps                           
t_f                                     = sampler_hold_fractions(sample.t_t,params                     );

%% Sample New Mean Molecule States
% [sample.mu_D,sample.mu_A,sample.mu_acc_rate_MH,sample.mu_acc_rate_HMC]=...
%                                               sampler_update_emission(sample.t_s,t_f,...
%                                               sample.mu_D,sample.mu_A,...
%                                               sample.mu_acc_rate_MH,...
%                                               sample.mu_acc_rate_HMC,...
%                                               params );
[sample.mu_D,sample.mu_A,...
 sample.mu_acc_rate_MH,...
 sample.mu_acc_rate_HMC,sample.acc_rate_flip]=...
                                              sampler_update_emission(sample.t_s,t_f,...
                                              sample.mu_D,sample.mu_A,...
                                              sample.mu_acc_rate_MH,...
                                              sample.mu_acc_rate_HMC,...
                                              sample.acc_rate_flip,...
                                              params );


% sample.mu_D         = params.ground.mu_D;
% sample.mu_A         = params.ground.mu_A;
% sample.mu_acc_rate  = [0 realmin];
%% Sample New Traj
[sample.t_t,sample.t_s]            = sampler_sample_traj(sample.t_t,sample.t_s,sample.mu_D,sample.mu_A,...
                                                         params.mu_back_D,params.mu_back_A,...
                                                         sample.IT,B_prob,t_f,params   );
%% Debugging rates and probs
%  sample.t_t                      = params.ground.t_t;
%  sample.t_s                      = params.ground.t_s;

%% Remove Virtual Jumps
[sample.t_t, sample.t_s]           = sampler_drop_jumps(sample.t_t,sample.t_s,params                   );

%% Calculate Holding Times After REMOVING the virtual jumps
t_d                                = sampler_hold_times(sample.t_t,params                              );

%% Transition Counts       
TCM                                = sampler_transit_count(sample.t_t,sample.t_s,params                );

%% Transition Probs
[sample.P,sample.IT]               = sampler_transit_prob(TCM,params                                   );

%% EscRate
sample.escrate                     = sampler_escrate(sample.t_s,TCM,t_d,params                         ); 

%% Rate Matrix
 sample.Q                          = sampler_full_rate(sample.P,sample.escrate,params                  );

t_f                                = sampler_hold_fractions(sample.t_t,params                     );
 
%% Calculate MAP
    
 sample.MAP =  calculate_MAP_est(sample.t_t,t_d,t_f,sample.t_s,...
                 sample.mu_D,params.mu_back_D,...
                 sample.mu_A,params.mu_back_A,...
                 sample.P,sample.IT,sample.Q,params,TCM);