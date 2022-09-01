function P_group = sampler_group_prob(B,IT,group_st,n_st ,flag)
%l_obs.
%k_group
%n_state
p_flag = flag;
switch p_flag
    case 0
    keyboard
    tag1 = (n_st(1,2)-n_st(1,1))+1;
    P_group = cell(1,length(group_st));
    switch tag1
        case 1    
        P_group{1} = IT';
        case 2
        P_group{1} = reshape(B',2^(tag1),1).*kron(IT',[1;1]);
        case 3
        E          = (kron(B,B).*repmat(kron(eye(2),[1 1]),2,1))  ;
        
        P_group{1} = [sum(E(1:2,:)) sum(E(3:4,:))]'.*kron(IT',[1;1;1;1]);
    end
                                  
    for k =2:length(group_st)
    
    tag = (n_st(k,2)-n_st(k,1))+1;

    switch tag
        case 1 
        keyboard
        P_group{k} = B;
        P_group{k} = repmat(B,2^(n_st(k-1,2)-n_st(k-1,1)),1);
        case 2
%         keyboard
        P          = kron(kron([1;1],eye(2)),[1 1]).*(multi_kron(B,n_st,k));
        P_group{k} = [sum(P(1:length(P(1,:))/2,:));sum(P((length(P(1,:))/2+1):end,:))];
        P_group{k} = repmat(P_group{k},2^((n_st(k-1,2)-n_st(k-1,1))),1);
        case 3
        keyboard
        P          = kron(kron([1;1],eye(4)),[1 1]).*(multi_kron(B,n_st,k));
        P_group{k} = [sum(P(1:length(P(1,:))/2,:));sum(P((length(P(1,:))/2+1):end,:))];
        P_group{k} = repmat(P_group{k},2^((n_st(k-1,2)-n_st(k-1,1))),1);
    end

    end
    case 1
        P_group = B;
end
tt = 4;
