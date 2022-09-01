function params = chainer_init_params(opts)

%% Units
params.units    =   opts.units   ;

%% Observation Parameters


params.Int_D      =   opts.Int_D     ;
params.Int_A      =   opts.Int_A     ;

params.t_left   =   opts.t_left  ;

params.t_right  =   opts.t_right ;

params.T_i      =   opts.T_i     ;

params.T_f      =   opts.T_f     ;


%% Hyperparameters
       
params.M        =   opts.M       ;
 
params.alpha    =   opts.alpha   ;
 
params.beta     =   opts.beta    ;

params.eta      =   opts.eta     ;   

params.psi_D    =   opts.psi_D   ;
params.phi_D    =   opts.phi_D   ;

params.psi_A    =   opts.psi_A   ;
params.phi_A    =   opts.phi_A   ;


% params.nu_D     =   opts.nu_D   ;
% params.chi_D    =   opts.chi_D   ;
% 
% params.nu_A     =   opts.nu_A   ;
% params.chi_A    =   opts.chi_A   ;


params.alpha_prop_D    =   opts.alpha_prop_D   ;
params.alpha_prop_A    =   opts.alpha_prop_A   ;

% params.theta_prop_D    =   opts.theta_prop_D   ;
% params.omega_prop_A    =   opts.omega_prop_A   ;

params.mu_back_D     = opts.mu_back_D;
params.mu_back_A     = opts.mu_back_A;

params.HMC_eps = opts.HMC_eps;
params.HMC_L = opts.HMC_L;
%% -----Debugging RATES and PROBS

if isfield(opts,'ground')
    
params.ground.t_t      = opts.ground.t_t;

params.ground.t_s      = opts.ground.t_s;

params.ground.Q        = opts.ground.Q;

params.ground.escrate  = opts.ground.escrate;

params.ground.mu_D     =opts.ground.mu_D;

params.ground.mu_A     = opts.ground.mu_A;
end

%% aux sampler's params

%  MCMC stride

params.i_skip   =   1            ;

