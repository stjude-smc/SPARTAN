function idx = chainer_alltimes_indices(t_t,params)

% l for observations
% k for panels
% n for times
% after additon of virtual jumps are added
% distribution of observation indices over all_times windows

idx         = cell(1,length(t_t));
obs_left      = params.t_left;
obs_right     = params.t_right;

for k = 1:length(t_t)
    if ~(abs(k-length(t_t))<eps)
        for l=1:length(params.obs)    
 
            %%% ----CASE2------
            if  t_t(k)<obs_left(l)&&t_t(k+1)<obs_right(l)&& t_t(k+1)>obs_left(l)
            
            idx{k} = [idx{k};l];
        
            %%% ----CASE3------
            elseif  t_t(k)>obs_left(l)&&t_t(k+1)>obs_right(l)&&t_t(k)<obs_right(l)
            idx{k} = [idx{k};l];

       
            %%% ----CASE4------
            elseif t_t(k)<obs_left(l) && t_t(k+1)>obs_right(l)
            idx{k} = [idx{k};l];

            %%% ----CASE5------
            elseif t_t(k)>obs_left(l)&&t_t(k+1)<obs_right(l)
            idx{k} = [idx{k};l];
        
            end

%             if isempty(idx{l})
%                  keyboard
%             end
        end 
    else
        for l=1:length(params.obs)    
            if  t_t(k)<obs_left(l) && params.T_f>obs_right(l)
            idx{k} = [idx{k};l];
        
        %%% ----CASE2------
            elseif  t_t(k)>obs_left(l)&& params.T_f>obs_right(l)&& t_t(k)<obs_right(l)
            idx{k} = [idx{k};l];

              
            end
%             if isempty(idx{l})
%                  keyboard
%             end
        end 
    end
end
 
% for k =1:length(t_t)
%     if isempty(idx{k})
%     	keyboard
%     end
% end
%  
 uu=0;

 

 
 
