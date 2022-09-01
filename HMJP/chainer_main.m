function chain = chainer_main(chain_init,d_length,opts,flag_status,flag_visual)
% to init:
% chain = chainer_main([]   ,  0, opts, true, []  );
% to expand:
% chain = chainer_main(chain,+25, []  , true, true);
% to reduce:
% chain = chainer_main(chain,-10, []  , true, []  );

rng('shuffle')


% init chain --------------------------------------------------------------
if d_length == 0

    % MCMC
    chain.params = chainer_init_params(opts);
    chain.length = 1;
    chain.sizeGB = nan;      % current chain memory size
    chain.sample = [];
    

    chain.sample = chainer_init_sample(chain.params,opts);   
    % history
    chain.i      = cast( chain.sample.i      , 'uint64') ;
    chain.mu_D   = cast( chain.sample.mu_D, 'single') ;
    chain.mu_A   = cast( chain.sample.mu_A, 'single') ;
%     chain.mu_back_D   = cast( chain.sample.mu_back_D, 'single') ;
%     chain.mu_back_A   = cast( chain.sample.mu_back_A, 'single') ;
    chain.t_t    = cast( chain.sample.t_t    , 'single') ;
    chain.t_s    = cast( chain.sample.t_s    , 'int16' ) ;
    chain.escrate= cast( chain.sample.escrate, 'single') ;
    chain.P      = cast( chain.sample.P      , 'single') ;
    chain.IT     = cast( chain.sample.IT     , 'single') ;
    chain.Q      = cast( chain.sample.Q      , 'single') ;
    chain.MAP    = cast( chain.sample.MAP     , 'single') ;

    
    if flag_status
        disp('CHAINER: chain initiated');
    end
    

% expand chain ------------------------------------------------------------
elseif d_length > 0

    chain.params = chain_init.params;
    chain.length = chain_init.length + d_length;
    chain.sizeGB = nan;
    chain.sample = chain_init.sample;
    chain.i       = [ chain_init.i      zeros( 1,d_length       , 'like', chain_init.i)];
    chain.mu_D    = [ chain_init.mu_D nan( chain.params.M, d_length,'like',chain_init.mu_D  )];
    chain.mu_A    = [ chain_init.mu_A nan( chain.params.M, d_length,'like',chain_init.mu_A  )];
%     chain.mu_back_D    = [ chain_init.mu_back_D nan( 1, d_length,'like',chain_init.mu_back_D  )];
%     chain.mu_back_A    = [ chain_init.mu_back_A nan( 1, d_length,'like',chain_init.mu_back_A  )];
    chain.t_t     = [ chain_init.t_t    ; cell( d_length, 1                             )]; 
    chain.t_s     = [ chain_init.t_s    ; cell( d_length, 1                             )]; 
    chain.escrate = [ chain_init.escrate; nan(  d_length,  chain.params.M ,'like',chain_init.escrate  )];
    chain.IT      = [ chain_init.IT     ; nan(  d_length,  chain.params.M ,'like',chain_init.IT  )];
    chain.P       = cat(3, chain_init.P ,  nan(   chain.params.M, chain.params.M,d_length )); 
    chain.Q       = cat(3, chain_init.Q ,  nan(   chain.params.M, chain.params.M,d_length )); 
    chain.MAP      = [ chain_init.MAP   ; nan(  d_length, 1,'like',chain_init.MAP)];

    
    if flag_visual
        Gim = chainer_visualize([],chain);
    end
    
    %---------------------------- expand chain
    r = chain_init.length+1;
    while r <= chain.length
        
        chain.sample = sampler_update(chain.sample,chain.params);
        
        if mod(chain.sample.i,chain.params.i_skip) == 0
            
            chain.i(r)           = chain.sample.i      ;
            chain.mu_D(:,r)      = chain.sample.mu_D   ;
            chain.mu_A(:,r)      = chain.sample.mu_A   ;
%             chain.mu_back_D(1,r) = chain.sample.mu_back_D   ;
%             chain.mu_back_A(1,r) = chain.sample.mu_back_A   ;
            chain.t_t{r}         = chain.sample.t_t    ;
            chain.t_s{r}         = chain.sample.t_s    ;
            chain.escrate(r,:)   = chain.sample.escrate;
            chain.IT(r,:)        = chain.sample.IT     ;
            chain.P(:,:,r)       = chain.sample.P      ;
            chain.Q(:,:,r)       = chain.sample.Q      ;
            chain.MAP(r,:)        = chain.sample.MAP     ;

            if flag_visual
                chainer_visualize(Gim,chain);
            end
            
            if flag_status
                disp([  'i = ', num2str(chain.sample.i,'%d'),... 
                    ' - MH acc = ' , num2str( chain.sample.mu_acc_rate_MH(1)/chain.sample.mu_acc_rate_MH(2) * 100 ,'%#6.2f') , ' %', ...
                    ' - HMC acc = ' , num2str( chain.sample.mu_acc_rate_HMC(1)/chain.sample.mu_acc_rate_HMC(2) * 100 ,'%#6.2f') , ' %',...        
                    ]);
            end
            
            r = r+1;
        end
    end    

    if flag_status
        
        disp('CHAINER: chain expanded');
        
    end


% reduce chain ------------------------------------------------------------
elseif d_length < 0

    d_length = -d_length;

    chain.params = chain_init.params;
    chain.length = d_length;
    chain.sizeGB = nan;
    chain.sample = chain_init.sample;
    
    ind = mod(chain_init.length,d_length)+(floor(chain_init.length/d_length)*(1:d_length));

    chain.i       = chain_init.i(ind)          ;
    chain.mu_D    = chain_init.mu_D(:,ind)  ;
    chain.mu_A    = chain_init.mu_D(:,ind)  ;
%     chain.mu_back_D = chain_init.mu_back_D(1,ind)  ;
%     chain.mu_back_A = chain_init.mu_back_A(1,ind)  ;
    for j =1:length(ind)
    chain.t_t{j}  = chain_init.t_t{ind(j)}        ;
    chain.t_s{j}     = chain_init.t_s{ind(j)}        ;
    end
    chain.escrate = chain_init.escrate(ind,:)  ;
    chain.IT      = chain_init.IT(ind,:)       ;
    chain.P       = chain_init.P(:,:,ind)      ;
    chain.Q       = chain_init.Q(:,:,ind)      ;
    chain.MAP     = chain_init.MAP(ind,:)       ;

    
    chain.params.i_skip = double(chain.i(2)-chain.i(1));
    
    
    if flag_status
        disp('CHAINER: chain reduced');
    end
end


%% book-keeping
chain.sizeGB = get_sizeGB(chain);               % mem size


end





%% auxiliary functions

function sizeGB = get_sizeGB(chain)
    sizeGB = whos( inputname(1) );
    sizeGB = sizeGB.bytes/1024^3;
end


