function P_multi = multi_kron(P,n_st,k)
% This function is responsible for kronocker product of the same matrix
% multiple times
PP1      =  P;
for j=1:(n_st(k,2)-n_st(k,1))
    PP1  = kron(P,PP1) ;
end
 P_multi = (PP1);
