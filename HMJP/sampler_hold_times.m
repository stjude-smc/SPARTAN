function hold_times = sampler_hold_times(t_t,params)

% l for observations
% k for panels
% n for times
% after additon of virtual jumps are added

hold_times  = [diff(t_t);params.T_f-(t_t(end))]; %column

% length(hold_times) = length(t_t) = length(t_s)
%% ----------------------DEBUGGING PURPOSES-------------------------------
if abs(sum(hold_times)-(params.T_f-params.T_i)) > 10*eps('single')
    
     keyboard
    
end


