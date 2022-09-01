function  TCM  = sampler_transit_count(t_t,t_s,params)



TCM                         =  zeros(params.M + 1,params.M);


TCM( params.M+1 ,t_s(1) )   = 1;

for n = 2 : length(t_s)
    TCM( t_s(n-1), t_s(n) ) = TCM( t_s(n-1), t_s(n) ) +1;
end





if ~all(diag(TCM)==0)
  keyboard
end