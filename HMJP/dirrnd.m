function p =dirrnd(rs)

p = randg(rs);
p = p./sum(p,2);
% p = randg(ab);
%p = bsxfun(@rdivide,p,sum(p,2));
 
 
% check for underflows and turn to stick-breakings if so
ind = isnan(p(:,1));
if any(ind)
    disp('Dir underflows')
     keyboard
    K = size(rs,2);
    p(ind,:) = zeros(sum(ind),K);
    for k=1:K
        p(ind,k) = betarnd(rs(ind,k),sum(rs(ind,k+1:K),2)+realmin).*(1-sum(p(ind,:),2));
    end
     
end
