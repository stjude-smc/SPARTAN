function [P,IT] = sampler_transit_prob(TCM,params)

%% Base Distribution

Diag                 = diag(ones(params.M,1));
logDiag              = logical( Diag(:) );
     vv              = ones(params.M);
     vv( logDiag )   = zeros(params.M,1);
     base            = (1/params.M)*(vv);  

     
%% Main Transition Probability Matrix Update
    P                = dirrnd( TCM(1:params.M,1:params.M) +base );

%% Initial Transition Probaility Matrix Update

    IT               = dirrnd(TCM(params.M+1,1)+(1/(params.M))*(ones(1,params.M)));

