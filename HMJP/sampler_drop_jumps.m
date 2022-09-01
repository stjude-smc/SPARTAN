function [t_t, t_s]  = sampler_drop_jumps(t_t,t_s,params)

drop_ind = [true;logical(diff(t_s))];
t_s      = t_s(drop_ind);
t_t      = t_t(drop_ind);

