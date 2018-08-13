function eQ = dtAdjustedQ( Q, td, classidx )
% Calculate rate matrix Q corrected for missed events with dead time td.
% classidx is a vector of class numbers for each state in Q.
% Qaa is the submatrix of Q involving only states in class a.
% Qab has transitions from states in class a to states in class b.
% Qaz has transitions from class a to states in *all other* classes.
% Qac has transitions from class a to states in classes other than a and b.
% See Qin 1996 and Milescu (2006) Biophys J (91), pg. 1156..
% See QtoQe() function in max_ll_util.cpp in QuB for implementation.

eQ = zeros(size(Q));
nClasses = max(classidx);

for aa=1:nClasses
    % Indices of states in current class for getting submatrices of Q.
    a  = classidx==aa;
    na = classidx~=aa;
    
    % Pre-compute factors used in every iteration
    Raa = ( eye(sum(na)) - expm(Q(na,na)*td) )  *  Q(na,na)^-1;  %eq. 18
    qra = Q(a,na) * Raa * Q(na,a);  
    
    for bb=1:nClasses
        b  = classidx==bb;
        nb = classidx~=bb;
        c  = na & nb;
        
        if a==b
            eQ(a,a) = Q(a,a) - qra;
        else
            eQ(a,b) = expm(td * qra) * ...
                  (  Q(a,b) - (Q(a,c) * ( eye(sum(c)) - expm(Q(c,c)*td) ) * ...
                     Q(c,c)^-1 * Q(c,b))  );
            %FIXME: handle two-class case, in which c is empty!
        end
    end %for each final state
end %for each starting state

% Renormalize so rows to sum zero.
% FIXME: should this not be needed if the code above is accurate???
I = logical(eye(size(eQ)));
eQ(I) = 0;
eQ(I) = -sum(eQ,2);

end %function dtAdjustedQ
