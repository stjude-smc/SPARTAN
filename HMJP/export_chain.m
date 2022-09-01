function export_chain(chain,fr, dr , nameit,frmt)


maxsample   = (length(chain.i));
minsample   = fix( fr*maxsample );
mem         = chain.sizeGB;
params      = chain.params;



if strcmp(frmt,'mat')

    for samplenum   = minsample:dr: maxsample

    t_s{samplenum-minsample+1,:}         = chain.t_s{samplenum,:};
    t_t{samplenum-minsample+1,:}         = chain.t_t{samplenum,:};
    IT(samplenum-minsample+1,:)          = chain.IT(samplenum,:);
    mu_D(:,samplenum-minsample+1)        = chain.mu_D(:,samplenum);
    mu_A(:,samplenum-minsample+1)        = chain.mu_A(:,samplenum);
    P(:,:,samplenum-minsample+1)         = chain.P(:,:,samplenum);
    Q(:,:,samplenum-minsample+1)         = chain.Q(:,:,samplenum);
    MAP(samplenum-minsample+1,1)         = chain.MAP(samplenum,:);
    end 

% cd '/home/zkilic/Dropbox (ASU)/zeliha_presselab/Project3_w_back_2_trace/data' ;
save([nameit,'.mat'],'mem','params','t_s','t_t','IT','mu_D','MAP','mu_A','P','Q','-v7.3');
%v7-3 is for data >2gb
disp([nameit,'.mat','--mission_accomplished']);
% cd '/home/zkilic/Dropbox (ASU)/zeliha_presselab/Project3_w_back_2_trace' ;
else
    disp('Unknown export format!')
end