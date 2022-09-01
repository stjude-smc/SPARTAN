function MAP = calculate_MAP_est(t_t,d,t_f,t_s,mu_D,mu_back_D,mu_A,mu_back_A,P,IT,Q,params,TCM)

%% Part1 - Contribution from the likelihood 

log_ld_s_D      = nan(length(t_t),1);
log_ld_s_A      = nan(length(t_t),1);

clear means_s_temp log_ld_s_k

    tStart = tic;
    dt = params.t_right(1)-params.t_left(1);

    mu_D_s_temp        = reshape(mu_D(t_s),1,1,length(t_s)); % obs*1*length(t_t)
    log_ld_s_D         = sum(-dt*(mu_back_D+sum(mu_D_s_temp.*t_f,3))+...
                               params.Int_D'.*log(dt*(mu_back_D+sum(mu_D_s_temp.*t_f,3)))-...
                               gammaln(params.Int_D'+1));   % 1 by M
                           
    mu_A_s_temp        = reshape(mu_A(t_s),1,1,length(t_s)); % obs*1*length(t_t)
    log_ld_s_A         = sum(-dt*(mu_back_A+sum(mu_A_s_temp.*t_f,3))+...
                               params.Int_A'.*log(dt*(mu_back_A+sum(mu_A_s_temp.*t_f,3)))-...
                               gammaln(params.Int_A'+1));   % 1 by M                       
 
    tElapsed            = toc(tStart);
    
%     tStart_2 = tic;
% 
%     for k=1:length(t_t)
%         
%     means_s_temp        = repmat(reshape(params.mean_st(t_s),1,1,length(t_s)),1,1,1);
%     log_ld_s_k_1        = -0.5*params.prec*sum((sum(means_s_temp.*t_f,3)-params.obs').^2);
% 
%     end
%     
%     tElapsed_2 = toc(tStart_2);


    MAP_1               = log_ld_s_D + log_ld_s_A;
        
%% Part2 - Contribution from the state transitions 
    
    logDiag             = logical(diag(ones(params.M,1)));
    P(logDiag )         = ones(params.M,1);
    MAP_2               = sum(sum(TCM(1:params.M,:).*log(P),2));
 


%% Part3 - Contribution from the holding time
    
    Q_d                = -diag(Q);

    MAP_3              = -sum( d./(Q_d(t_s) + log( Q_d(t_s) ) ) );
    
%% Part4 - Contribution from the initial state transition

    MAP_4              = sum( TCM(end,:).*log(IT) );
    
%% Part5 - Contribution from the embeded chain matrix

    MAP_5              = sum( sum( (1-1)*log(P) ) );

%% Part6 - Contribution from the initial transition probbaility matrix

    MAP_6              = sum( sum( (params.alpha/2-1)*log(IT) ) );

%% Part7 - Contribution coming from teh escape rates

    MAP_7              = sum( (params.eta/2-1)*log(Q_d) - params.eta*params.beta*Q_d/2 );
    
%% Part8 - Contribution coming from the emission rates-DONOR

    MAP_8              = sum( (params.phi_D/2-1)*log(mu_D) - params.phi_D*params.psi_D*mu_D/2 );
    
%% Part9 - Contribution coming from the BACKGROUND emission rates-DONOR

    MAP_9              = 0;%sum( (params.chi_D/2-1)*log(mu_back_D) - params.chi_D*params.nu_D*mu_back_D/2 );
    
%% Part8 - Contribution coming from the emission rates-ACCEPTOR

    MAP_10             = sum( (params.phi_A/2-1)*log(mu_A) - params.phi_A*params.psi_A*mu_A/2 );
    
%% Part9 - Contribution coming from the BACKGROUND emission rates-ACCEPTOR

    MAP_11             =0;% sum( (params.chi_A/2-1)*log(mu_back_A) - params.chi_A*params.nu_A*mu_back_A/2 );
    
    
%% FINAL MAP
MAP = MAP_1+ MAP_2+ MAP_3+ MAP_4+ MAP_5+ MAP_6+ MAP_7+ MAP_8+ MAP_9+ MAP_10+ MAP_11;


