save_file = [mfilename,'.mat'];
format compact
%% simulate data
tag              = 1;
tag_case         = 2;
% [observation,units,ground] = generate_synthetic_data(tag,true);
load('simulated_data_1_alpha_1.mat')
% observation        = load('FRET_data_tag_1_dr_0.011_bin_1000000000.mat');
% 
% observation.t_left  = [observation.t_edges_final(1:end-1) observation.t_edges_final(end)-(1e+8)]*10^-9;
% observation.t_right = [observation.t_edges_final(1)+(1e+8) observation.t_edges_final(2:end)]*10^-9;
% units.time = 's';
% units.space= 'photons';


%% init chain

opts.Int_D       = observation.Int_D;
opts.Int_A       = observation.Int_A;

opts.t_left      = observation.t_left;

opts.t_right     = observation.t_right;

opts.units       = units;

opts.M           = 3;

opts.alpha       = 2;
 
opts.beta        = 0.01;

opts.eta         = 4;

opts.phi_D         = 1; % mean of the mean_st DONOR
opts.psi_D         = sum(opts.Int_D)/(ground.T_f-ground.T_i);%2000; %prec for mean_st  DONOR

opts.phi_A         = 1; % mean of the mean_st ACCEPTOR
opts.psi_A         = sum(opts.Int_A)/(ground.T_f-ground.T_i);% 2000; %prec for mean_st  ACCEPTOR

% opts.nu_D          =  sum(opts.Int_D)/(ground.T_f-ground.T_i);%2000; % mean of the mean_st DONOR
% opts.chi_D         = 1; %prec for mean_st  DONOR
% 
% opts.nu_A          = sum(opts.Int_A)/(ground.T_f-ground.T_i);% 2000; % mean of the mean_st ACCEPTOR
% opts.chi_A         = 1; %prec for mean_st  ACCEPTOR

opts.alpha_prop_D  = sum(opts.Int_D)/(ground.T_f-ground.T_i);% 2000;
opts.alpha_prop_A  = sum(opts.Int_A)/(ground.T_f-ground.T_i);% 2000;

opts.theta_prop_D  = sum(opts.Int_D)/(ground.T_f-ground.T_i);%2000;
opts.omega_prop_A  = sum(opts.Int_A)/(ground.T_f-ground.T_i);%2000;


opts.mu_back_D   = ground.mu_back_D;
opts.mu_back_A   = ground.mu_back_A;

opts.HMC_eps = 0.01;
opts.HMC_L   = 50;

%% ---------DEBUGGING RATES and PROBS
if ~isempty(ground)
opts.ground.t_t         = ground.traj_t;
opts.ground.t_s         = ground.traj_s;
opts.ground.Q           = ground.Q;
opts.ground.escrate     = -diag(ground.Q)';
opts.ground.mu_D        = ground.mu_D;
opts.ground.mu_A        = ground.mu_A;
end


opts.T_i                = ground.T_i;
opts.T_f                = ground.T_f;
    
chain_length            = 20;

chain = chainer_main([],0,opts,true,[]);

clear opts


%% expand chain

while true
                     
    chain = chainer_main(chain,chain_length,[],true,true);

    if chain.sizeGB > 1.0
    
    chain = chainer_main(chain,-fix(chain.length/2),[],true,[]);
    
    end
    save(save_file)
    disp(['SAVED: ', save_file])
    
    plot_escrate;
    
    plot_diag_prob;
    
    plot_trajectories;
    
    plot_state_level;
drawnow
end

