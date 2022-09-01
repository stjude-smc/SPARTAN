function chain = reanalyze_samples_true(name)


%% load the data
%eg name = 'tag_all_test_traj_known_10_17601';
cd  '/home/zkilic/Dropbox (ASU)/zeliha_presselab/Project3_w_meanst_poisson_2_trace_w_BACK_EXP/data' ; 

data = load( [name '.mat'] );

cd  '/home/zkilic/Dropbox (ASU)/zeliha_presselab/Project3_w_meanst_poisson_2_trace_w_BACK_EXP' ; 

s_data           = size(data.IT);
chain.length     = s_data(1);
mu_D               = nan(2,1);
mu_A               = nan(2,1);

chain.mu_D       = nan(2,chain.length  );
chain.mu_back_D  = nan(1,chain.length  );
chain.mu_A       = nan(2,chain.length  );
chain.mu_back_A  = nan(1,chain.length  );
chain.t_t        = cell(chain.length,1 );
chain.t_s        = cell(chain.length,1 );
chain.Q          = nan(2,2,chain.length);
chain.IT         = nan(chain.length,2  );
for j=1:chain.length
    
    mu_D  =  data.mu_D(:,j);
    mu_A  =  data.mu_A(:,j);
    
    qu  = -diag(data.Q(:,:,j));
    
    if ((data.mu_D(1,j)>data.mu_D(2,j))|(data.mu_A(1,j)>data.mu_A(2,j)))       
    chain.mu_D(1,j)              = mu_D(2);
    chain.mu_D(2,j)              = mu_D(1);
    chain.mu_A(1,j)              = mu_A(2);
    chain.mu_A(2,j)              = mu_A(1);
          escrate(1,j)           = qu(2);
          escrate(2,j)           = qu(1);
    idx                          = data.t_s{j}==1;
    idx_p                        = data.t_s{j}==2;

    data.t_s{j}(idx,1)           = 2;
    data.t_s{j}(idx_p,1)         = 1;
    chain.t_t{j}                 = data.t_t{j};
    chain.t_s{j}                 = data.t_s{j};  
    chain.Q(1,2,j)               = escrate(1,j);
    chain.Q(1,1,j)               = -escrate(1,j);
    chain.Q(2,1,j)               = escrate(2,j);
    chain.Q(2,2,j)               = -escrate(2,j);
    chain.IT(j,1 )               = data.IT(j,2 );      
    chain.IT(j,2 )               = data.IT(j,1 );      
    chain.mu_back_D(:,j)         = data.mu_back_D(:,j);
    chain.mu_back_A(:,j)         = data.mu_back_A(:,j);
    
    else
    chain.t_t{j}                 = data.t_t{j};
    chain.t_s{j}                 = data.t_s{j}; 
    chain.mu_D(:,j)              = data.mu_D(:,j);
    chain.mu_A(:,j)              = data.mu_A(:,j);
    chain.mu_back_D(:,j)         = data.mu_back_D(:,j);
    chain.mu_back_A(:,j)         = data.mu_back_A(:,j);
    chain.Q(1,1,j)               = data.Q(1,1,j);
    chain.Q(1,2,j)               = data.Q(1,2,j);
    chain.Q(2,1,j)               = data.Q(2,1,j);
    chain.Q(2,2,j)               = data.Q(2,2,j);
    chain.IT(j,:)                = data.IT(j,: );
    end     

end
    chain.P = data.P;
    chain.sizeGB   = data.mem;
    chain.params   = data.params;
   

chain.MAP  = nan(length(chain.t_t),1);

for i=1:length(data.t_t)

    t_f        = sampler_hold_fractions(chain.t_t{i},chain.params                    );
     
    TCM        = sampler_transit_count(chain.t_t{i},chain.t_s{i},chain.params       );  
    
    d          = sampler_hold_times(chain.t_t{i},chain.params                       );
    
  chain.MAP(i) =  calculate_MAP_est(chain.t_t{i},d,t_f,data.t_s{i},...
                 chain.mu_D(:,i),chain.mu_back_D(:,i),...
                 chain.mu_A(:,i),chain.mu_back_A(:,i),...
                 chain.P(:,:,i),chain.IT(i,:),chain.Q(:,:,i),chain.params,TCM);
end

cd '/home/zkilic/Dropbox (ASU)/zeliha_presselab/Project3_w_meanst_poisson_2_trace_w_BACK_EXP/data' ;
save([name,'_newsq_','.mat'],'chain','-v7.3');
%v7-3 is for data >2gb
disp([name,'.mat','--mission_accomplished']);
cd '/home/zkilic/Dropbox (ASU)/zeliha_presselab/Project3_w_meanst_poisson_2_trace_w_BACK_EXP' ;
