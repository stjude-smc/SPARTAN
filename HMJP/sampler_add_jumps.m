function [all_t, all_s] = sampler_add_jumps(t_t,t_s,Q,omega,params)
% k :for panels

%% Join Trajectory and Obs ----> Create Panels


all_times                = t_t;

all_states               = t_s;

% Never EVER!!! No holding state jump times can go beyond T_f!!!!!!!!!!!!!!


 
 all_durations          = [diff(all_times); params.T_f-all_times(end)];
 loc_event              = [];
    for k = 1 : length(all_states)       
        enu_event  = poissrnd((omega + Q(all_states(k),all_states(k)))*all_durations(k));       
        loc_event  = [loc_event;all_times(k) + all_durations(k)*rand(enu_event,1)]; %column       
    end
    
    temp_times = sort([loc_event;all_times],'ascend');
    
    all_s      = add_pts(all_times,all_states,temp_times,params); 

    all_t      = temp_times;
         
        

    

