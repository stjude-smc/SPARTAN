function[B_prob,omega] = sampler_uniformize(Q_rate,opts)
%% Uniformizing the Continuous Time Process

omega         =  (opts.alpha)*max(abs(diag(Q_rate)));

dim           =  size(Q_rate);

B_prob        =  eye(dim(1)) + (Q_rate/(omega));

