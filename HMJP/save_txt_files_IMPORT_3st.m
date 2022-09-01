%% script for saving the files in .txt format by Zeliha Kilic 

 
 %% M=3
%%% load files 100% duty cycle, generate_synthetic_data_multi_duty_cycle(1,true)

load('simulated_data_1_alpha_1.mat')
Int_D = observation.Int_D;
Int_A = observation.Int_A;
mu_back = [ground.mu_back_D,ground.mu_back_A];

 save('example_donor_trace_100percent_duty_cycle_M_3.txt','Int_D','-ascii');
 save('example_acceptor_trace_100percent_duty_cycle_M_3.txt','Int_A','-ascii');
 save('example_background_100percent_duty_cycle_M_3.txt','mu_back','-ascii');

