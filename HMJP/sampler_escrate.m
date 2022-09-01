function escrate = sampler_escrate(t_s,TCM,t_d,params)


sum_TCM   = sum(TCM,1)';

% Total holding times at visited states

total_t   = nan(params.M,1);

for j=1:params.M
    total_t(j) = sum( t_d(abs(t_s-j)<eps) );
end
% 
% total_t11   = [ sum( t_d(abs(t_s-1)<eps) );    
%                 sum( t_d(abs(t_s-2)<eps) )];

          
escrate   = gamrnd(sum_TCM + params.eta/2, 1./( total_t + (params.eta .* params.beta)*0.5 ) )';

if ~isrow(escrate)
    keyboard
end






% sum_TCM   = sum(TCM,1)';
% 
% % Total holding times at visited states
% 
% total_t   = [ sum( t_d(abs(t_s-1)<eps) );    
%               sum( t_d(abs(t_s-2)<eps) )];
% 
% escrate   = gamrnd(sum_TCM + params.eta/2, 1./( total_t + (params.eta .* params.beta)*0.5 ) )';
% 
% if ~isrow(escrate)
%     keyboard
% end
% 
% 

