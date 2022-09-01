function data_mjp_tag_13 = define_params(name)

%% load data
load([name,'.mat']);



params.obs         = observation.obs;

params.t_left      = observation.t_left;

params.t_right     = observation.t_right;

params.prec        = observation.prec;

params.mean_st     = observation.mean_st;

params.units       = units;

params.M           = length(params.mean_st);

params.alpha       = 2;
 
params.beta        = 0.03;

params.eta         = 4;


params.ground.t_t         = ground.traj_t;
params.ground.t_s         = ground.traj_s;
params.ground.Q           = ground.Q;
params.ground.escrate     = -diag(ground.Q)';
params.T_i                = ground.T_i;
params.T_f                = ground.T_f;

data_mjp_tag_13.params    = params;

