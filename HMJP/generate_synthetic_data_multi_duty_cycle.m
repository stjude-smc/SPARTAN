function [observation,units,ground] = generate_synthetic_data_multi_duty_cycle(tag,show_demo)

addpath('fig')

if nargin<1
    show_demo = true;   
    tag       = 1;
end

%% set up-Units

% units
 units.time    = 's';
 units.space   = 'nm';
 units.Int       = 'photons';
 

 switch tag
     case 1    
      dt_1       = 0.05;%100 ms window
     alpha     = 1;
      tau_f    = .1;
     case 2
       dt_1      = 0.05;%1000 ms window
     alpha     = 1;
      tau_f    = .05;
     case 3
      dt_1       = 0.05;%100 ms window
     alpha     = 1;
      tau_f    = .8;
     case 4
       dt_1      =0.05;%1000 ms window
     alpha     = 1;
      tau_f    = .5;  
     case 5
       dt_1      = 0.05;%1000 ms window
     alpha     = 1;
      tau_f    = 1/7;
 end
 
%% Observation Window Determination


obs_wind_left  = (0.1 :dt_1 :4-(0.08))'; %same size as obs(COLUMN)
obs_wind_right = obs_wind_left + alpha*dt_1; %0 .01 is the exposure (COLUMN)
            

%% set up-parameters 
  %decide  here observatoion windows 
 %observation windows first then T_i(<<<<T_L1 of the obseravtion) is small then  and 
 %T_f>>>> Last side of the observation Window 

       
Q            =  (1/tau_f)*[-2.2 2.2;
                           2.1 -2.1]; 
s_init        =  1; 

T_i           =  obs_wind_left(1)      - (0.05); %
                      
T_f           =  obs_wind_right(end)   +  0.1;  % 


% 0.1 is critical this is critical, when you choose small again
% again nan appears after interp1

mu_D            =  [800;2000];%photons/sec  mean state(1 2 3) = (1 2 5)(units of space)
mu_A            =  [300;700];%photons/sec  mean state(1 2 3) = (1 2 5)(units of space)

mu_back_D       =  200;
mu_back_A       =  200;


%% simulate data and such

M             =  zeros(2,2);

for i=1:2
    for j=1:2
        if ~isequal(i,j)
        M(i,j) = Q(i,j)/(-Q(i,i));
        else
        M(i,i) = 0;
        end
    end
end

     time          = T_i;

     stID          = s_init;
     
     trace_t       = time;
     
     trace_s       = stID;


while (time < T_f)
          
    t_d         = log(rand) / (Q(stID,stID));
    
    stID        = sample_categorical(M(stID,:));%it samples based onthe rows of the transition matrix
    
    time        = time + t_d;

    trace_t     = [trace_t; time]; 
    
    trace_s     = [trace_s; stID];
    
end
 
trace_t_full    = trace_t;

trace_s_full    = trace_s;


% last_ind        = find(trace_t>T_f,1);

trace_t         = trace_t(1:end-1);
trace_s         = trace_s(1:end-1);



%% ground-truth

 ground.Q           = Q;
 ground.traj_t      = trace_t;
 ground.traj_s      = trace_s;
 ground.T_i         = T_i;
 ground.T_f         = T_f;
 ground.mu_D        = mu_D;
 ground.mu_A        = mu_A;
 ground.mu_back_D   = mu_back_D; %photons/sec
 ground.mu_back_A   = mu_back_A; %photons/sec

%% observation extraction
% spit out the observations and the observation windows

%%% All observations

[tx,sx]    =  stairs(trace_t_full,trace_s_full);

%% Fixing the holding states: (Jump Times, Holding States)

sleft     =  nan(length(trace_s_full(:,1)),1);
tleft     =  nan(length(trace_t_full(:,1)),1);%column

for n     =  1:length(trace_s_full(:,1))
    
sleft(n)  =  sx(2*n-1); % Holding State IDS
tleft(n)  =  tx(2*n-1); % Jump times to the holding states
    
end


sright    = nan(length(trace_t_full(:,1))-1,1);%column
tright    = nan(length(trace_t_full(:,1))-1,1);%column

for n=1:length(trace_t_full(:,1))-1
    
sright(n) = sx(2*n);
tright(n) = tx(2*n);
    
end

  
%% Holding States



Int_D               = nan(1,length(obs_wind_right)); % Holding State Preallocation(ROW)
Int_A               = nan(1,length(obs_wind_right)); % Holding State Preallocation(ROW)

options.T_f         = T_f;
options.T_i         = T_i;


%% Finding the indexes contributing to an observation
 
        

        all_times                = sort([obs_wind_left;obs_wind_right;trace_t],'ascend');
        
        all_states               = add_pts(trace_t,trace_s,all_times,options);
                   
        all_t_d                  = [diff(all_times);nan];
        
        ff                       = nan(length(obs_wind_left),1);
        
        idx                      = cell(1,length(obs_wind_left));
        
        for n   = 1:length(obs_wind_left)
            
        idx{n}  = find((all_times>=obs_wind_left(n)) & all_times<obs_wind_right(n));
        
        if ((abs(sum(all_t_d(idx{n}))-(obs_wind_right(1)-obs_wind_left(1))))>10*eps) %recall 0.01 is
                                                                 %the window size
            keyboard
        end
        
        ff(n)   = sum(all_t_d(idx{n})) ;
        
        end 
                
%-------------------------------------------------------------------------- 


%% Average Calculation  during integration time



for j =1:length(obs_wind_right)
    Int_D(j) = poissrnd(mu_back_D*(obs_wind_right(1)-obs_wind_left(1))+(sum(all_t_d(idx{j}).*mu_D(all_states(idx{j})))));
    Int_A(j) = poissrnd(mu_back_A*(obs_wind_right(1)-obs_wind_left(1))+(sum(all_t_d(idx{j}).*mu_A(all_states(idx{j})))));
end

Int_D = max(Int_D,realmin);
Int_A = max(Int_A,realmin);

observation.Int_D   = Int_D;
observation.Int_A   = Int_A;
observation.t_left  = obs_wind_left';
observation.t_right = obs_wind_right';

save(['simulated_tag_less_precise_' num2str(tag),'_dt_',num2str(dt_1),'_alpha_',num2str(alpha),'.mat'],'ground','observation','units','-v7.3')
if isequal(length(tleft),2)

        disp('no transition occured!');
else
    %%
if show_demo

    
    % ---------------------------------------------------------------------
    fig = figure(1);
    fig.Name = 'Gillespie Trajectory';
    clf
  axes1 =   subplot(2,1,1);
  axes2 =   subplot(2,1,2);

%% Gillespie Trajectory
axes(axes1)

p1(1)      =  plot(tx(1:end-1),sx(1:end-1),'-mo',...
  'LineWidth',1);

hold on

p1(2)      =  plot(tleft(1:end-1)',sleft(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

p1(3)      =  plot(tright(1:end-1)',sright(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.9 1 .9],...
    'MarkerSize',5);

p1(3)      =  plot(tright(end)',sright(end)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

p1(1).LineWidth = 1 ;
line(tright(1:end-1)'.*[1;1],[0,4],'linestyle',':','color','k');
line(T_i.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);


xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.05*(T_f-T_i);
xlim_end =  T_f + 0.05*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','\sigma_{3}','',''})

title('Gillespie Trajectory')


hold off


%% Points from the trajectory
axes(axes2)

pp(1)      =  plot(trace_t',trace_s','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

line(T_i.*[1 1],[0,6],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,6],'linestyle','--','color','k','linewidth',2);
line(trace_t(1:end-1).*[1 1],[0,6],'linestyle',':','color','k');


pp(1).LineWidth = 1 ;

xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','\sigma_{3}','',''})

title('Holding States and Jump Times from Gillespie')

hold off
% exportgraphics data_set_1_gillespie_traj -pdf
%% Observations 
    fig = figure(2);
    fig.Name = 'Observations from the Gillespie trajectory';
    clf
  axes1 =   subplot(2,1,1);
  axes2 =   subplot(2,1,2);
  
%%% ONLY Gillespie
%%% Trajectory-----------------------------------------------------------------
  
axes(axes1)

ppp(1)      =  plot(tx(1:end-1),sx(1:end-1),'-mo');

ppp(1).LineWidth = 1 ;

hold on

ppp(2)      =  plot(tleft(1:end-1)',sleft(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

ppp(3)      =  plot(tright(1:end-1)',sright(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.9 1 .9],...
    'MarkerSize',5);

ppp(4)      =  plot(tright(end)',sright(end)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

line(tright(1:end-1)'.*[1;1],[0,4],'linestyle',':','color','k');
line(T_i.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);


xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','','',''})

ylabel({'States'})

title('Gillespie Trajecory')



%------------------------------------------------------------------------------
%%% ONLY Observations 

axes(axes2)
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(obs_wind_right)
    
l_A =  line([obs_wind_left(j) obs_wind_right(j)],[Int_A(j) Int_A(j)],'color','r','LineWidth',5);
l_D =  line([obs_wind_left(j) obs_wind_right(j)],[Int_D(j) Int_D(j)],'color','g','LineWidth',5);

hold on

end


hold on


xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;


xlim([xlim_0 xlim_end])




l1 = line((obs_wind_left)'.*[1 ;1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle',':','color','g','LineWidth',1);

l2 = line(obs_wind_right'.*[1; 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle',':','color','b','LineWidth',1);

% l3 = line([T_i,T_f],[mean_state(1);(mean_state(2))]'.*[1;1],'linestyle','-','color','c','LineWidth',1);

line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);

line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);

legend([l_D(1) l_A(1) l1(1) l2(1)],'Donor channel','Acceptor channel','{observation window left end point}','{observation window right end point}')

hold off

xlabel('Time (s)')
ylabel({'Intensity','(photon)'})
ylim([min([Int_A Int_D])-10,max([Int_A Int_D])+10])


xlim([xlim_0 xlim_end])

title('Simulated Observations')
box on
% export_fig data_set_1_gillespie_traj_obs -pdf

%% Combined demo

    fig = figure(3);
    fig.Name = 'State of the molcule and intensities';
    clf
  
    col_S = [0 0 1];   % main sytem color
    col_D = [0 1 0];   % main donor color
    col_A = [1 0 0];   % main acceptor color
    col_m = [0 1 1];   % aux molecule color


    M = length(mu_D);
    tn_bnd           =  [T_i;obs_wind_left(2:end-1);T_f];% [t]

    
    
    
    subplot(4,1,1)
    h = stairs(tx(1:end-1),sx(1:end-1));
%     hold on 
%     plot_stairs(tn_bnd,sn);
    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_S;
    ylabel({'State of','molecule'})
    xlim([xlim_0 xlim_end])
    ylim([0 M+1])
    set(gca,'YTick',1:M,'YTickLabel',[repmat('\sigma_',M,1),num2str((1:M)')])
    legend('State of the molecule','Extrap. state of the molecule','location','North','orientation','horizontal')
    legend boxoff
    box off
   
    
    subplot(4,1,[2 3])
    h(1) = stairs(tn_bnd,Int_D);
    hold on
    h(2) = stairs(tn_bnd,Int_A);

    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_D;
    h(2).Color = col_A;
    ylabel({'Intensities',['(',units.Int,')']})
    xlim([xlim_0 xlim_end])
    legend('Donor channel','Acceptor channel','location','North','orientation','horizontal')
    legend boxoff
    box off
    
    subplot(4,1,4)
    h(1) = stairs(tn_bnd,100*Int_A./(Int_D+Int_A));
    hold on
    h(2) = stairs(trace_t,(100*mu_A(trace_s)./(mu_D(trace_s)+mu_A(trace_s)))');
    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_S;
    h(2).Color = col_m;
    ylabel({'FRET','efficiency (%)'})
    xlabel(['Time (',units.time,')'])
    ylim([-20 75])
    xlim([xlim_0 xlim_end])
    legend('Apparent','True','location','North','orientation','horizontal')
    legend boxoff
    box off
  
% exportgraphics data_set_1_analyzed_obs -pdf
  
  


end
end
%% Part 2

%% Observation Window Determination
alpha = 0.75*alpha;

obs_wind_left  = (0.1 :dt_1 :4-(0.08))'; %same size as obs(COLUMN)
obs_wind_right = obs_wind_left + alpha*dt_1; %0 .01 is the exposure (COLUMN)
            
T_i           =  obs_wind_left(1)      - (0.05); %
                      
T_f           =  obs_wind_right(end)   +  0.1;  % 

%% ground-truth

 ground.Q           = Q;
 ground.traj_t      = trace_t;
 ground.traj_s      = trace_s;
 ground.T_i         = T_i;
 ground.T_f         = T_f;
 ground.mu_D        = mu_D;
 ground.mu_A        = mu_A;
 ground.mu_back_D   = mu_back_D; %photons/sec
 ground.mu_back_A   = mu_back_A; %photons/sec
 
%% Holding States



Int_D               = nan(1,length(obs_wind_right)); % Holding State Preallocation(ROW)
Int_A               = nan(1,length(obs_wind_right)); % Holding State Preallocation(ROW)

options.T_f         = T_f;
options.T_i         = T_i;


%% Finding the indexes contributing to an observation
 
        

        all_times                = sort([obs_wind_left;obs_wind_right;trace_t],'ascend');
        
        all_states               = add_pts(trace_t,trace_s,all_times,options);
                   
        all_t_d                  = [diff(all_times);nan];
        
        ff                       = nan(length(obs_wind_left),1);
        
        idx                      = cell(1,length(obs_wind_left));
        
        for n   = 1:length(obs_wind_left)
            
        idx{n}  = find((all_times>=obs_wind_left(n)) & all_times<obs_wind_right(n));
        
        if ((abs(sum(all_t_d(idx{n}))-(obs_wind_right(1)-obs_wind_left(1))))>10*eps) %recall 0.01 is
                                                                 %the window size
            keyboard
        end
        
        ff(n)   = sum(all_t_d(idx{n})) ;
        
        end 
                
%-------------------------------------------------------------------------- 


%% Average Calculation  during integration time



for j =1:length(obs_wind_right)
    Int_D(j) = poissrnd(mu_back_D*(obs_wind_right(1)-obs_wind_left(1))+(sum(all_t_d(idx{j}).*mu_D(all_states(idx{j})))));
    Int_A(j) = poissrnd(mu_back_A*(obs_wind_right(1)-obs_wind_left(1))+(sum(all_t_d(idx{j}).*mu_A(all_states(idx{j})))));
end

Int_D = max(Int_D,realmin);
Int_A = max(Int_A,realmin);

observation.Int_D   = Int_D;
observation.Int_A   = Int_A;
observation.t_left  = obs_wind_left';
observation.t_right = obs_wind_right';

save(['simulated_tag_less_precise_' num2str(tag),'_dt_',num2str(dt_1),'_alpha_',num2str(alpha),'.mat'],'ground','observation','units','-v7.3')

if isequal(length(tleft),2)

        disp('no transition occured!');
else
    %%
if show_demo

    
    % ---------------------------------------------------------------------
    fig = figure(4);
    fig.Name = 'Gillespie Trajectory';
    clf
  axes1 =   subplot(2,1,1);
  axes2 =   subplot(2,1,2);

%% Gillespie Trajectory
axes(axes1)

p1(1)      =  plot(tx(1:end-1),sx(1:end-1),'-mo',...
  'LineWidth',1);

hold on

p1(2)      =  plot(tleft(1:end-1)',sleft(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

p1(3)      =  plot(tright(1:end-1)',sright(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.9 1 .9],...
    'MarkerSize',5);

p1(3)      =  plot(tright(end)',sright(end)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

p1(1).LineWidth = 1 ;
line(tright(1:end-1)'.*[1;1],[0,4],'linestyle',':','color','k');
line(T_i.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);


xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.05*(T_f-T_i);
xlim_end =  T_f + 0.05*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','\sigma_{3}','',''})

title('Gillespie Trajectory')


hold off


%% Points from the trajectory
axes(axes2)

pp(1)      =  plot(trace_t',trace_s','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

line(T_i.*[1 1],[0,6],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,6],'linestyle','--','color','k','linewidth',2);
line(trace_t(1:end-1).*[1 1],[0,6],'linestyle',':','color','k');


pp(1).LineWidth = 1 ;

xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','\sigma_{3}','',''})

title('Holding States and Jump Times from Gillespie')

hold off
% export_fig data_set_2_gillespie_traj -pdf

%% Observations 
    fig = figure(5);
    fig.Name = 'Observations from the Gillespie trajectory';
    clf
  axes1 =   subplot(2,1,1);
  axes2 =   subplot(2,1,2);
  
%%% ONLY Gillespie
%%% Trajectory-----------------------------------------------------------------
  
axes(axes1)

ppp(1)      =  plot(tx(1:end-1),sx(1:end-1),'-mo');

ppp(1).LineWidth = 1 ;

hold on

ppp(2)      =  plot(tleft(1:end-1)',sleft(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

ppp(3)      =  plot(tright(1:end-1)',sright(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.9 1 .9],...
    'MarkerSize',5);

ppp(4)      =  plot(tright(end)',sright(end)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

line(tright(1:end-1)'.*[1;1],[0,4],'linestyle',':','color','k');
line(T_i.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);


xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','','',''})

ylabel({'States'})

title('Gillespie Trajecory')



%------------------------------------------------------------------------------
%%% ONLY Observations 

axes(axes2)
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(obs_wind_right)
    
l_A =  line([obs_wind_left(j) obs_wind_right(j)],[Int_A(j) Int_A(j)],'color','r','LineWidth',5);
l_D =  line([obs_wind_left(j) obs_wind_right(j)],[Int_D(j) Int_D(j)],'color','g','LineWidth',5);

hold on

end


hold on


xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;


xlim([xlim_0 xlim_end])




l1 = line((obs_wind_left)'.*[1 ;1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle',':','color','g','LineWidth',1);

l2 = line(obs_wind_right'.*[1; 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle',':','color','b','LineWidth',1);

% l3 = line([T_i,T_f],[mean_state(1);(mean_state(2))]'.*[1;1],'linestyle','-','color','c','LineWidth',1);

line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);

line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);

legend([l_D(1) l_A(1) l1(1) l2(1)],'Donor channel','Acceptor channel','{observation window left end point}','{observation window right end point}')

hold off

xlabel('Time (s)')
ylabel({'Intensity','(photon)'})
ylim([min([Int_A Int_D])-10,max([Int_A Int_D])+10])


xlim([xlim_0 xlim_end])

title('Simulated Observations')
box on
% exportgraphics data_set_2_gillespie_traj_obs -pdf

%% Combined demo

    fig = figure(6);
    fig.Name = 'State of the molcule and intensities';
    clf
  
    col_S = [0 0 1];   % main sytem color
    col_D = [0 1 0];   % main donor color
    col_A = [1 0 0];   % main acceptor color
    col_m = [0 1 1];   % aux molecule color


    M = length(mu_D);
    tn_bnd           =  [T_i;obs_wind_left(2:end-1);T_f];% [t]

    
    
    
    subplot(4,1,1)
    h = stairs(tx(1:end-1),sx(1:end-1));
%     hold on 
%     plot_stairs(tn_bnd,sn);
    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_S;
    ylabel({'State of','molecule'})
    xlim([xlim_0 xlim_end])
    ylim([0 M+1])
    set(gca,'YTick',1:M,'YTickLabel',[repmat('\sigma_',M,1),num2str((1:M)')])
    legend('State of the molecule','Extrap. state of the molecule','location','North','orientation','horizontal')
    legend boxoff
    box off
   
    
    subplot(4,1,[2 3])
    h(1) = stairs(tn_bnd,Int_D);
    hold on
    h(2) = stairs(tn_bnd,Int_A);

    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_D;
    h(2).Color = col_A;
    ylabel({'Intensities',['(',units.Int,')']})
    xlim([xlim_0 xlim_end])
    legend('Donor channel','Acceptor channel','location','North','orientation','horizontal')
    legend boxoff
    box off
    
    subplot(4,1,4)
    h(1) = stairs(tn_bnd,100*Int_A./(Int_D+Int_A));
    hold on
    h(2) = stairs(trace_t,(100*mu_A(trace_s)./(mu_D(trace_s)+mu_A(trace_s)))');
    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_S;
    h(2).Color = col_m;
    ylabel({'FRET','efficiency (%)'})
    xlabel(['Time (',units.time,')'])
    ylim([-20 75])
    xlim([xlim_0 xlim_end])
    legend('Apparent','True','location','North','orientation','horizontal')
    legend boxoff
    box off
  
%   exportgraphics data_set_2_analyzed_obs -pdf

  


end
end 
 
%% Part 3
%% Observation Window Determination
alpha = (4/3)*(1/2)*alpha;

obs_wind_left  = (0.1 :dt_1 :4-(0.08))'; %same size as obs(COLUMN)
obs_wind_right = obs_wind_left + alpha*dt_1; %0 .01 is the exposure (COLUMN)
            
T_i           =  obs_wind_left(1)      - (0.05); %
                      
T_f           =  obs_wind_right(end)   +  0.1;  % 

%% ground-truth

 ground.Q           = Q;
 ground.traj_t      = trace_t;
 ground.traj_s      = trace_s;
 ground.T_i         = T_i;
 ground.T_f         = T_f;
 ground.mu_D        = mu_D;
 ground.mu_A        = mu_A;
 ground.mu_back_D   = mu_back_D; %photons/sec
 ground.mu_back_A   = mu_back_A; %photons/sec
 
%% Holding States



Int_D               = nan(1,length(obs_wind_right)); % Holding State Preallocation(ROW)
Int_A               = nan(1,length(obs_wind_right)); % Holding State Preallocation(ROW)

options.T_f         = T_f;
options.T_i         = T_i;


%% Finding the indexes contributing to an observation
 
        

        all_times                = sort([obs_wind_left;obs_wind_right;trace_t],'ascend');
        
        all_states               = add_pts(trace_t,trace_s,all_times,options);
                   
        all_t_d                  = [diff(all_times);nan];
        
        ff                       = nan(length(obs_wind_left),1);
        
        idx                      = cell(1,length(obs_wind_left));
        
        for n   = 1:length(obs_wind_left)
            
        idx{n}  = find((all_times>=obs_wind_left(n)) & all_times<obs_wind_right(n));
        
        if ((abs(sum(all_t_d(idx{n}))-(obs_wind_right(1)-obs_wind_left(1))))>10*eps) %recall 0.01 is
                                                                 %the window size
            keyboard
        end
        
        ff(n)   = sum(all_t_d(idx{n})) ;
        
        end 
                
%-------------------------------------------------------------------------- 


%% Average Calculation  during integration time



for j =1:length(obs_wind_right)
    Int_D(j) = poissrnd(mu_back_D*(obs_wind_right(1)-obs_wind_left(1))+(sum(all_t_d(idx{j}).*mu_D(all_states(idx{j})))));
    Int_A(j) = poissrnd(mu_back_A*(obs_wind_right(1)-obs_wind_left(1))+(sum(all_t_d(idx{j}).*mu_A(all_states(idx{j})))));
end

Int_D = max(Int_D,realmin);
Int_A = max(Int_A,realmin);

observation.Int_D   = Int_D;
observation.Int_A   = Int_A;
observation.t_left  = obs_wind_left';
observation.t_right = obs_wind_right';

save(['simulated_tag_less_precise_' num2str(tag),'_dt_',num2str(dt_1),'_alpha_',num2str(alpha),'.mat'],'ground','observation','units','-v7.3')

 
if isequal(length(tleft),2)

        disp('no transition occured!');
else
    %%
if show_demo

    
    % ---------------------------------------------------------------------
    fig = figure(7);
    fig.Name = 'Gillespie Trajectory';
    clf
  axes1 =   subplot(2,1,1);
  axes2 =   subplot(2,1,2);

%% Gillespie Trajectory
axes(axes1)

p1(1)      =  plot(tx(1:end-1),sx(1:end-1),'-mo',...
  'LineWidth',1);

hold on

p1(2)      =  plot(tleft(1:end-1)',sleft(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

p1(3)      =  plot(tright(1:end-1)',sright(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.9 1 .9],...
    'MarkerSize',5);

p1(3)      =  plot(tright(end)',sright(end)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

p1(1).LineWidth = 1 ;
line(tright(1:end-1)'.*[1;1],[0,4],'linestyle',':','color','k');
line(T_i.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);


xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.05*(T_f-T_i);
xlim_end =  T_f + 0.05*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','\sigma_{3}','',''})

title('Gillespie Trajectory')


hold off


%% Points from the trajectory
axes(axes2)

pp(1)      =  plot(trace_t',trace_s','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

line(T_i.*[1 1],[0,6],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,6],'linestyle','--','color','k','linewidth',2);
line(trace_t(1:end-1).*[1 1],[0,6],'linestyle',':','color','k');


pp(1).LineWidth = 1 ;

xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','\sigma_{3}','',''})

title('Holding States and Jump Times from Gillespie')

hold off
% export_fig data_set_3_gillespie_traj -pdf

%% Observations 
    fig = figure(8);
    fig.Name = 'Observations from the Gillespie trajectory';
    clf
  axes1 =   subplot(2,1,1);
  axes2 =   subplot(2,1,2);
  
%%% ONLY Gillespie
%%% Trajectory-----------------------------------------------------------------
  
axes(axes1)

ppp(1)      =  plot(tx(1:end-1),sx(1:end-1),'-mo');

ppp(1).LineWidth = 1 ;

hold on

ppp(2)      =  plot(tleft(1:end-1)',sleft(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

ppp(3)      =  plot(tright(1:end-1)',sright(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.9 1 .9],...
    'MarkerSize',5);

ppp(4)      =  plot(tright(end)',sright(end)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

line(tright(1:end-1)'.*[1;1],[0,4],'linestyle',':','color','k');
line(T_i.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);


xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','','',''})

ylabel({'States'})

title('Gillespie Trajecory')



%------------------------------------------------------------------------------
%%% ONLY Observations 

axes(axes2)
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(obs_wind_right)
    
l_A =  line([obs_wind_left(j) obs_wind_right(j)],[Int_A(j) Int_A(j)],'color','r','LineWidth',5);
l_D =  line([obs_wind_left(j) obs_wind_right(j)],[Int_D(j) Int_D(j)],'color','g','LineWidth',5);

hold on

end


hold on


xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;


xlim([xlim_0 xlim_end])




l1 = line((obs_wind_left)'.*[1 ;1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle',':','color','g','LineWidth',1);

l2 = line(obs_wind_right'.*[1; 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle',':','color','b','LineWidth',1);

% l3 = line([T_i,T_f],[mean_state(1);(mean_state(2))]'.*[1;1],'linestyle','-','color','c','LineWidth',1);

line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);

line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);

legend([l_D(1) l_A(1) l1(1) l2(1)],'Donor channel','Acceptor channel','{observation window left end point}','{observation window right end point}')

hold off

xlabel('Time (s)')
ylabel({'Intensity','(photon)'})
ylim([min([Int_A Int_D])-10,max([Int_A Int_D])+10])


xlim([xlim_0 xlim_end])

title('Simulated Observations')
box on

% exportgraphics data_set_3_gillespie_traj_obs -pdf

%% Combined demo

    fig = figure(9);
    fig.Name = 'State of the molcule and intensities';
    clf
  
    col_S = [0 0 1];   % main sytem color
    col_D = [0 1 0];   % main donor color
    col_A = [1 0 0];   % main acceptor color
    col_m = [0 1 1];   % aux molecule color


    M = length(mu_D);
    tn_bnd           =  [T_i;obs_wind_left(2:end-1);T_f];% [t]

    
    
    
    subplot(4,1,1)
    h = stairs(tx(1:end-1),sx(1:end-1));
%     hold on 
%     plot_stairs(tn_bnd,sn);
    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_S;
    ylabel({'State of','molecule'})
    xlim([xlim_0 xlim_end])
    ylim([0 M+1])
    set(gca,'YTick',1:M,'YTickLabel',[repmat('\sigma_',M,1),num2str((1:M)')])
    legend('State of the molecule','Extrap. state of the molecule','location','North','orientation','horizontal')
    legend boxoff
    box off
   
    
    subplot(4,1,[2 3])
    h(1) = stairs(tn_bnd,Int_D);
    hold on
    h(2) = stairs(tn_bnd,Int_A);

    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_D;
    h(2).Color = col_A;
    ylabel({'Intensities',['(',units.Int,')']})
    xlim([xlim_0 xlim_end])
    legend('Donor channel','Acceptor channel','location','North','orientation','horizontal')
    legend boxoff
    box off
    
    subplot(4,1,4)
    h(1) = stairs(tn_bnd,100*Int_A./(Int_D+Int_A));
    hold on
    h(2) = stairs(trace_t,(100*mu_A(trace_s)./(mu_D(trace_s)+mu_A(trace_s)))');
    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_S;
    h(2).Color = col_m;
    ylabel({'FRET','efficiency (%)'})
    xlabel(['Time (',units.time,')'])
    ylim([-20 75])
    xlim([xlim_0 xlim_end])
    legend('Apparent','True','location','North','orientation','horizontal')
    legend boxoff
    box off
  
  
%   export_fig data_set_3_analyzed_obs -pdf
end
end
%% Part 4
%% Observation Window Determination
alpha = (1/2)*alpha;

obs_wind_left  = (0.1 :dt_1 :4-(0.08))'; %same size as obs(COLUMN)
obs_wind_right = obs_wind_left + alpha*dt_1; %0 .01 is the exposure (COLUMN)
            
T_i           =  obs_wind_left(1)      - (0.05); %
                      
T_f           =  obs_wind_right(end)   +  0.1;  % 

%% ground-truth

 ground.Q           = Q;
 ground.traj_t      = trace_t;
 ground.traj_s      = trace_s;
 ground.T_i         = T_i;
 ground.T_f         = T_f;
 ground.mu_D        = mu_D;
 ground.mu_A        = mu_A;
 ground.mu_back_D   = mu_back_D; %photons/sec
 ground.mu_back_A   = mu_back_A; %photons/sec
 
%% Holding States



Int_D               = nan(1,length(obs_wind_right)); % Holding State Preallocation(ROW)
Int_A               = nan(1,length(obs_wind_right)); % Holding State Preallocation(ROW)

options.T_f         = T_f;
options.T_i         = T_i;


%% Finding the indexes contributing to an observation
 
        

        all_times                = sort([obs_wind_left;obs_wind_right;trace_t],'ascend');
        
        all_states               = add_pts(trace_t,trace_s,all_times,options);
                   
        all_t_d                  = [diff(all_times);nan];
        
        ff                       = nan(length(obs_wind_left),1);
        
        idx                      = cell(1,length(obs_wind_left));
        
        for n   = 1:length(obs_wind_left)
            
        idx{n}  = find((all_times>=obs_wind_left(n)) & all_times<obs_wind_right(n));
        
        if ((abs(sum(all_t_d(idx{n}))-(obs_wind_right(1)-obs_wind_left(1))))>10*eps) %recall 0.01 is
                                                                 %the window size
            keyboard
        end
        
        ff(n)   = sum(all_t_d(idx{n})) ;
        
        end 
                
%-------------------------------------------------------------------------- 


%% Average Calculation  during integration time



for j =1:length(obs_wind_right)
    Int_D(j) = poissrnd(mu_back_D*(obs_wind_right(1)-obs_wind_left(1))+(sum(all_t_d(idx{j}).*mu_D(all_states(idx{j})))));
    Int_A(j) = poissrnd(mu_back_A*(obs_wind_right(1)-obs_wind_left(1))+(sum(all_t_d(idx{j}).*mu_A(all_states(idx{j})))));
end

Int_D = max(Int_D,realmin);
Int_A = max(Int_A,realmin);

observation.Int_D   = Int_D;
observation.Int_A   = Int_A;
observation.t_left  = obs_wind_left';
observation.t_right = obs_wind_right';

save(['simulated_tag_less_precise_' num2str(tag),'_dt_',num2str(dt_1),'_alpha_',num2str(alpha),'.mat'],'ground','observation','units','-v7.3')

 
if isequal(length(tleft),2)

        disp('no transition occured!');
else
    %%
if show_demo

    
    % ---------------------------------------------------------------------
    fig = figure(7);
    fig.Name = 'Gillespie Trajectory';
    clf
  axes1 =   subplot(2,1,1);
  axes2 =   subplot(2,1,2);

%% Gillespie Trajectory
axes(axes1)

p1(1)      =  plot(tx(1:end-1),sx(1:end-1),'-mo',...
  'LineWidth',1);

hold on

p1(2)      =  plot(tleft(1:end-1)',sleft(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

p1(3)      =  plot(tright(1:end-1)',sright(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.9 1 .9],...
    'MarkerSize',5);

p1(3)      =  plot(tright(end)',sright(end)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

p1(1).LineWidth = 1 ;
line(tright(1:end-1)'.*[1;1],[0,4],'linestyle',':','color','k');
line(T_i.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);


xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.05*(T_f-T_i);
xlim_end =  T_f + 0.05*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','\sigma_{3}','',''})

title('Gillespie Trajectory')


hold off


%% Points from the trajectory
axes(axes2)

pp(1)      =  plot(trace_t',trace_s','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

line(T_i.*[1 1],[0,6],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,6],'linestyle','--','color','k','linewidth',2);
line(trace_t(1:end-1).*[1 1],[0,6],'linestyle',':','color','k');


pp(1).LineWidth = 1 ;

xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','\sigma_{3}','',''})

title('Holding States and Jump Times from Gillespie')

hold off
% exportgraphics data_set_4_gillespie_traj -pdf

%% Observations 
    fig = figure(8);
    fig.Name = 'Observations from the Gillespie trajectory';
    clf
  axes1 =   subplot(2,1,1);
  axes2 =   subplot(2,1,2);
  
%%% ONLY Gillespie
%%% Trajectory-----------------------------------------------------------------
  
axes(axes1)

ppp(1)      =  plot(tx(1:end-1),sx(1:end-1),'-mo');

ppp(1).LineWidth = 1 ;

hold on

ppp(2)      =  plot(tleft(1:end-1)',sleft(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

ppp(3)      =  plot(tright(1:end-1)',sright(1:end-1)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.9 1 .9],...
    'MarkerSize',5);

ppp(4)      =  plot(tright(end)',sright(end)','o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.3 .3 .3],...
    'MarkerSize',5);

line(tright(1:end-1)'.*[1;1],[0,4],'linestyle',':','color','k');
line(T_i.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);
line(T_f.*[1 1],[0,4],'linestyle','--','color','k','linewidth',2);


xlabel('Time (s)')
ylabel('States')
ylim([0,4])

xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;

xlim([xlim_0 xlim_end])
yticks([1 2 3 4 ])
yticklabels({'\sigma_{1}','\sigma_{2}','','',''})

ylabel({'States'})

title('Gillespie Trajecory')



%------------------------------------------------------------------------------
%%% ONLY Observations 

axes(axes2)
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(obs_wind_right)
    
l_A =  line([obs_wind_left(j) obs_wind_right(j)],[Int_A(j) Int_A(j)],'color','r','LineWidth',5);
l_D =  line([obs_wind_left(j) obs_wind_right(j)],[Int_D(j) Int_D(j)],'color','g','LineWidth',5);

hold on

end


hold on


xlim_0   =  T_i - 0.07*(T_f-T_i);
xlim_end =  T_f + 0.07*(T_f-T_i) ;


xlim([xlim_0 xlim_end])




l1 = line((obs_wind_left)'.*[1 ;1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle',':','color','g','LineWidth',1);

l2 = line(obs_wind_right'.*[1; 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle',':','color','b','LineWidth',1);

% l3 = line([T_i,T_f],[mean_state(1);(mean_state(2))]'.*[1;1],'linestyle','-','color','c','LineWidth',1);

line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);

line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);

legend([l_D(1) l_A(1) l1(1) l2(1)],'Donor channel','Acceptor channel','{observation window left end point}','{observation window right end point}')

hold off

xlabel('Time (s)')
ylabel({'Intensity','(photon)'})
ylim([min([Int_A Int_D])-10,max([Int_A Int_D])+10])


xlim([xlim_0 xlim_end])

title('Simulated Observations')
box on

% exportgraphics data_set_4_gillespie_traj_obs -pdf

%% Combined demo

    fig = figure(9);
    fig.Name = 'State of the molcule and intensities';
    clf
  
    col_S = [0 0 1];   % main sytem color
    col_D = [0 1 0];   % main donor color
    col_A = [1 0 0];   % main acceptor color
    col_m = [0 1 1];   % aux molecule color


    M = length(mu_D);
    tn_bnd           =  [T_i;obs_wind_left(2:end-1);T_f];% [t]

    
    
    
    subplot(4,1,1)
    h = stairs(tx(1:end-1),sx(1:end-1));
%     hold on 
%     plot_stairs(tn_bnd,sn);
    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_S;
    ylabel({'State of','molecule'})
    xlim([xlim_0 xlim_end])
    ylim([0 M+1])
    set(gca,'YTick',1:M,'YTickLabel',[repmat('\sigma_',M,1),num2str((1:M)')])
    legend('State of the molecule','Extrap. state of the molecule','location','North','orientation','horizontal')
    legend boxoff
    box off
   
    
    subplot(4,1,[2 3])
    h(1) = stairs(tn_bnd,Int_D);
    hold on
    h(2) = stairs(tn_bnd,Int_A);

    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_D;
    h(2).Color = col_A;
    ylabel({'Intensities',['(',units.Int,')']})
    xlim([xlim_0 xlim_end])
    legend('Donor channel','Acceptor channel','location','North','orientation','horizontal')
    legend boxoff
    box off
    
    subplot(4,1,4)
    h(1) = stairs(tn_bnd,100*Int_A./(Int_D+Int_A));
    hold on
    h(2) = stairs(trace_t,(100*mu_A(trace_s)./(mu_D(trace_s)+mu_A(trace_s)))');
    line(T_i.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    line(T_f.*[1 1],[min([Int_A Int_D])-10,max([Int_A Int_D])+10],'linestyle','--','color','k','linewidth',2);
    h(1).Color = col_S;
    h(2).Color = col_m;
    ylabel({'FRET','efficiency (%)'})
    xlabel(['Time (',units.time,')'])
    ylim([-20 75])
    xlim([xlim_0 xlim_end])
    legend('Apparent','True','location','North','orientation','horizontal')
    legend boxoff
    box off
  
  
%   exportgraphics data_set_4_analyzed_obs -pdf



end
end




        