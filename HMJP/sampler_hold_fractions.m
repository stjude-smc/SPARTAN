function hold_fractions = sampler_hold_fractions(t_t,params)

% l for observations
% k for panels
% n for times
% after additon of virtual jumps are added

hold_times    = zeros( length(params.Int_D),1,length(t_t) );
% diff_hold_t   = [diff(t_t);params.T_f-params.T_i];
%%%%??

for k=1:length(t_t)
    if ~(abs(k-length(t_t))<eps)
        for l = 1:length(params.t_left)
            %%% ----CASE1------
            if     t_t(k+1)<=params.t_left(l)
            hold_times(l,1,k) = hold_times(l,1,k)+0;
        
            %%% ----CASE2------
            elseif  t_t(k)<=params.t_left(l)&&t_t(k+1)<params.t_right(l)&& t_t(k+1)>params.t_left(l)
            hold_times(l,1,k) =  hold_times(l,1,k) + (t_t(k+1)-params.t_left(l));
        
            %%% ----CASE3------
            elseif  t_t(k)>=params.t_left(l)&&t_t(k+1)>params.t_right(l)&&t_t(k)<params.t_right(l)
            hold_times(l,1,k) =  hold_times(l,1,k) + (params.t_right(l)-t_t(k));
            
            %%% ----CASE4------
            elseif t_t(k)<=params.t_left(l) && t_t(k+1)>params.t_right(l)
            hold_times(l,1,k) =  hold_times(l,1,k) + (params.t_right(l)-params.t_left(l));

            %%% ----CASE5------
            elseif t_t(k)>=params.t_left(l)&&t_t(k+1)<params.t_right(l)
            hold_times(l,1,k) = hold_times(l,1,k) + (t_t(k+1)-t_t(k));
        
            %%% ----CASE6------
            elseif t_t(k)>=params.t_right(l)    
            hold_times(l,1,k) = hold_times(l,1,k)+0;
        
            end

        end
    else

        for l = 1:length(params.t_left)
                
        %%% ----CASE1------
            if  t_t(k)<=params.t_left(l) && params.T_f>=params.t_right(l)
            hold_times(l,1,k) =  hold_times(l,1,k) + (params.t_right(l)-params.t_left(l));
        
        %%% ----CASE2------
            elseif  t_t(k)>=params.t_left(l)&& params.T_f>=params.t_right(l)&& t_t(k)<params.t_right(l)
            hold_times(l,1,k) =  hold_times(l,1,k) + (params.t_right(l)-t_t(k));
      
        %%% ----CASE3------
            elseif t_t(k)>=params.t_right(l)    
            hold_times(l,1,k) = hold_times(l,1,k)+0;
        
            end
            
            
        end
    end    
        
        
   
end

hold_fractions = hold_times./(params.t_right-params.t_left)';
hold_frac_reshaped = reshape(hold_fractions,length(params.Int_D),length(t_t));

if any((abs(sum(hold_frac_reshaped,2)-1)>eps('single')) )
    keyboard
end

if any(reshape(abs(sum(hold_fractions,3)-1)>eps('single'),1,[] ) )
    keyboard
end
uu = 9;

