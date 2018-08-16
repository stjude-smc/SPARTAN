function eQ = dtAdjustedQ( Q, td, classidx )
% Calculate rate matrix Q corrected for missed events with dead time td.
% classidx is a vector of class numbers for each state in Q.
% Qaa is the submatrix of Q involving only states in class a.
% Qab has transitions from states in class a to states in class b.
% Qaz has transitions from class a to states in *all other* classes.
% Qac has transitions from class a to states in classes other than a and b.
% See Qin 1996 and Milescu (2006) Biophys J (91), pg. 1156..
% See QtoQe() function in max_ll_util.cpp in QuB for implementation.

eQ = Q; %zeros(size(Q));
nClasses = max(classidx);

for aa=1:nClasses
    % Indices of states in current class for getting submatrices of Q.
    a  = classidx==aa;
    nota = classidx~=aa;
    
    % Pre-compute common factors for all end states
    % NOTE: eq. 18 uses Qaa and 16-17 use Qnota,nota. Which is it??
    Raa = ( eye(sum(nota)) - expm(Q(nota,nota)*td) )  *  Q(nota,nota)^-1;
    qra = Q(a,nota) * Raa * Q(nota,a);
    
    % Calculate eQaa.
    % Accounts for excursions out of class a that were too short to detect.
    eQ(a,a) = Q(a,a) - qra;
    
    % Calculate eQab for each target class b.
    % Accounts for indirect excursions a->c->a, without detecting dwell in c.
    for bb=1:nClasses
        if aa==bb, continue; end
        
        b = classidx==bb;
        notb = classidx~=bb;
        c = nota & notb;
        
        if sum(c)~=0
            escIndirect = Q(a,c) * ( eye(sum(c)) - expm(Q(c,c)*td) ) * Q(c,c)^-1 * Q(c,b);
        else
            % Can't escape indirectly if no third class exists.
            % This implementation matches QuB. Is it correct??
            escIndirect = zeros( sum(a), sum(b) );
        end
        eQ(a,b) = expm(td * qra) * ( Q(a,b) - escIndirect );
        
    end %for each final state
    
end %for each starting state


% Evidently, we dont want to renormalize eQ. If we do, eQ ends up being
% essentially the same as Q. Why??


end %function dtAdjustedQ
