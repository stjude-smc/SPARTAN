function sample = chainer_init_sample(params,opts)


%% Counter

sample.i = 0;

%% Transition probailtiy matrix initialization

%--------------------------------------------------------------------------

       Diag             = diag(ones(params.M,1));
    logDiag             = logical( Diag(:) );
     vv                 = ones(opts.M);
     vv( logDiag )      = zeros(params.M,1);
     base               = (1/params.M)*(vv);
      sample.P          = dirrnd(base);
          
         sample.IT      = dirrnd((1/(params.M))*ones(1,params.M));
         
%--------------------------------------------------------------------------         
%% MAP initialization         
%--------------------------------------------------------------------------
        sample.MAP      = 0;
  
%--------------------------------------------------------------------------         
%% Escape rate initialization         
%--------------------------------------------------------------------------

   sample.escrate       = gamrnd(params.eta/2,2./(opts.eta .* params.beta),[1 params.M]);
%--------------------------------------------------------------------------

         sample.Q       = -diag(sample.escrate) +  (sample.P).*(sample.escrate)';
%--------------------------------------------------------------------------

    
%% Initial State of the trajectory

%--------------------------------------------------------------------------

              s_i       = sample_categorical(sample.IT);
              
%--------------------------------------------------------------------------         
%% Mol state level initialization
%--------------------------------------------------------------------------
sample.mu_D          = (params.psi_D/params.phi_D)*randg(params.phi_D,[params.M,1]);

sample.mu_D          = sort(sample.mu_D);

sample.mu_A          = (params.psi_A/params.phi_A)*randg(params.phi_A,[params.M,1]);

sample.mu_A          = sort(sample.mu_A);

   sample.mu_acc_rate_MH    = [0 realmin];
   sample.mu_acc_rate_HMC   = [0 realmin];
   sample.acc_rate_flip     = [0 realmin];





%    sample.mu_D          = (params.psi_D/params.phi_D)*randg(params.phi_D,[params.M,1]);
%    sample.mu_A          = (params.psi_A/params.phi_A)*randg(params.phi_A,[params.M,1]);
%    sample.mu_acc_rate_MH    = [0 realmin];
%    sample.mu_acc_rate_HMC   = [0 realmin];

%-------------------------  -------------------------------------------------  

%% Trajectory initialization      

%--------------------------------------------------------------------------          
    
 [sample.t_t,sample.t_s] = chain_initial_sample_my(sample.Q,s_i,opts,false);

%             sample.t_t = params.ground.t_t; %ground truth change
%             sample.t_s = params.ground.t_s;



%% INITIALIZATION of TRAJECTORY
%%%%%%%%%%%%%%%%%%%%%%%%%%Initialization of the%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%      trajectory     %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[trace_t, trace_s]= chain_initial_sample_my(Q,s_i,opts,show_demo)


if nargin<1 
    
    [observation,~,~]     = generate_synthetic_data;
    
            Q             = (1/tau_f)*[-.05 .05;
                                        .1 -.1];                       
    s_i                   = 1;
    
    t_right               = observation.t_right;
    
    t_left                = observation.t_left;
    
    opts.T_i              = t_left(1)            ; %  
    
    opts.T_i              = t_right(end)         ;  % 
    
    show_demo             = true;
    
end


%--------------------------------------------------------------------------


        T_i     = opts.T_i;
        
        T_f     = opts.T_f;

%--------------------------------------------------------------------------

M               = zeros(opts.M,opts.M);

for i = 1:opts.M
    for j = 1:opts.M
        if ~isequal(i,j)
        M(i,j)  = Q(i,j)/(-Q(i,i));
        else
        M(i,i)  = 0;
        end
    end
end

     time       = T_i;

     stID       = s_i;
     
     trace_t    = time;
     
     trace_s    = stID;


while (time < T_f)
          
    holdTime    = log(rand) / (Q(stID,stID));
    
    stID        = sample_categorical(M(stID,:));% it samples based onthe 
                                                % rows of the transition
                                                % matrix  
    time        = time + holdTime;

    trace_t     = [trace_t; time]; 
    
    trace_s     = [trace_s; stID];
    
end

% last_ind        = find(trace_t>T_f,1);

trace_t         = trace_t(1:end-1);
trace_s         = trace_s(1:end-1);

%--------------------------------------------------------------------------
if show_demo

    [tx,sx]     =  stairs(trace_t,trace_s);

%------------PLOTTING Purposes---------------------------------------
    sleft       =  nan(length(trace_s(:,1)),1);
    
    tleft       =  nan(length(trace_t(:,1)),1);%column

    for n       =  1:length(trace_s(:,1))
    
    sleft(n)    =  sx(2*n-1); % Holding State IDS
   
    tleft(n)    =  tx(2*n-1); % Jump times to the holding states
    
    end


    sright      = nan(length(trace_t(:,1))-1,1);%column

    tright      = nan(length(trace_t(:,1))-1,1);%column

    for n=1:length(trace_t(:,1))-1
    
      sright(n) = sx(2*n);
      tright(n) = tx(2*n);
    
    end


%--------------------------------------------------------------------------


    
% -------------------------------------------------------------------------
    fig = figure(3);
    
    fig.Name = 'Initial Trajectory';
    clf
    axes1       = subplot(1,1,1);

    
    
axes(axes1)

p1(1)           =  plot(tx,sx,'-mo',...
  'LineWidth',1);

hold on

p1(2)           =  plot(tleft(1:end-1)',sleft(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);


p1(3)           =  plot(tright(1:end-1)',sright(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.9 1 .9],...
    'MarkerSize',5);

p1(4)           =  plot(tright(end)',sright(end)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

p1(1).LineWidth = 1 ;

line(tright(1:end-1).*[1 1],[0,4],'linestyle',':','color','k');
line(T_i.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);


xlabel('Time (s)')

ylabel('States')

ylim([0,4])

xlim_0          =  T_i - 0.05*(T_f-T_i);
xlim_end        =  T_f + 0.05*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','','',''})

title('Gillespie Trajectory')

end





