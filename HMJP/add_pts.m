function all_states = add_pts(t_times,t_states,all_times,params)

if length(t_times)>1
    
all_states = interp1(t_times,t_states,all_times,'previous','extrap');

else
    
t_times  = [t_times;(params.T_f+params.T_i)*0.5];

t_states = [t_states;t_states]; 

all_states = interp1(t_times,t_states,all_times,'previous','extrap');

end
    



