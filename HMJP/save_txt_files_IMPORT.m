%% script for saving the files in .txt format by Zeliha Kilic 

%% M=2
%%% load files 100% duty cycle, generate_synthetic_data_multi_duty_cycle(1,true)
% % 
% % load('simulated_tag_less_precise_1_dt_0.05_alpha_1.mat')
% % Int_D = observation.Int_D;
% % Int_A = observation.Int_A;
% % mu_back = [ground.mu_back_D,ground.mu_back_A];
% % 
% %  save('example_donor_trace_100percent_duty_cycle.txt','Int_D','-ascii');
% %  save('example_acceptor_trace_100percent_duty_cycle.txt','Int_A','-ascii');
% %  save('example_background_100percent_duty_cycle.txt','mu_back','-ascii');
% % 
% % clear Int_D Int_A mu_back
% % 
% % %%% load files 50% duty cycle
% % 
% % load('simulated_tag_less_precise_1_dt_0.05_alpha_0.5.mat')
% % Int_D = observation.Int_D;
% % Int_A = observation.Int_A;
% % mu_back = [ground.mu_back_D,ground.mu_back_A];
% % 
% %  save('example_donor_trace_50percent_duty_cycle.txt','Int_D','-ascii');
% %  save('example_acceptor_trace_50percent_duty_cycle.txt','Int_A','-ascii');
% %  save('example_background_50percent_duty_cycle.txt','mu_back','-ascii');
 
 
 
 %% M=3
%%% load files 100% duty cycle, generate_synthetic_data_multi_duty_cycle(1,true)

load('simulated_data_1_alpha_1.mat')
Int_D = observation.Int_D;
Int_A = observation.Int_A;
mu_back = [ground.mu_back_D,ground.mu_back_A];

 save('example_donor_trace_100percent_duty_cycle.txt','Int_D','-ascii');
 save('example_acceptor_trace_100percent_duty_cycle.txt','Int_A','-ascii');
 save('example_background_100percent_duty_cycle.txt','mu_back','-ascii');

clear Int_D Int_A mu_back

%%% load files 50% duty cycle

load('simulated_tag_less_precise_1_dt_0.05_alpha_0.5.mat')
Int_D = observation.Int_D;
Int_A = observation.Int_A;
mu_back = [ground.mu_back_D,ground.mu_back_A];

 save('example_donor_trace_50percent_duty_cycle.txt','Int_D','-ascii');
 save('example_acceptor_trace_50percent_duty_cycle.txt','Int_A','-ascii');
 save('example_background_50percent_duty_cycle.txt','mu_back','-ascii');