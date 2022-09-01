function [t_t,t_s] = chainer_FFBS(t_t,t_s,mean_st,t_f,params,P_gr,group_st,n_st,l_obs)

%l observ.
%k groups
%n states

for k = 1:length(group_st)
    TT            = cell(1,n_st(k,2)-n_st(k,1)+1);
    states        = [1 2];
        for kk = 1:n_st(k,2)-n_st(k,1)+1
         TT{kk}= states;
        end
    st_id{k}   = allcomb(TT{:});     
end

for k = 1:length(group_st)   
means_s_temp = reshape(mean_st(st_id{k}),1,length(st_id{k}),length(t_s(n_st(k,1):n_st(k,2)))); % obs*1*length(t_t)
L_part{k}    = -0.5*params.prec*sum((sum(means_s_temp.*t_f(l_obs(k,1):l_obs(k,2),1,n_st(k,1):n_st(k,2)),3)-params.obs(l_obs(k,1):l_obs(k,2))').^2);
L_part{k}    = L_part{k} -max(L_part{k});
L_part{k}    = exp( L_part{k} )';
L_part{k}    = max(L_part{k},realmin);
end


%% Filter forward
    L_part{1} = L_part{1} .* P_gr{1};
    L_part{1} = L_part{1} / sum(L_part{1});
for k = 2 : +1 : length(group_st)
    L_part{k} = L_part{k} .* ( P_gr{k}' * L_part{k-1} );
    L_part{k} = L_part{k} / sum(L_part{k});
end

%% Sample backward

s_id{length(group_st)}  = sample_categorical(L_part{length(group_st)});
t_s_gr{length(group_st)}= st_id{length(group_st)}(s_id{length(group_st)},:)';
for k=length(group_st)-1:-1:1
    s_id{k}             = sample_categorical( L_part{k}.*log(P_gr{k+1}));
    t_s_gr{k}           = st_id{k}(s_id{k},:)';
end

%% Convert super states
t_s = [];
for k = 1:length(group_st)
    t_s = [t_s;t_s_gr{k}];
end









