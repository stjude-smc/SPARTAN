function strans=transmat_uncor

%--- cut of traces after being in the 3rd state for e.g. 10000ms 
MAX_THREE_DWELL=10000; 

ntraces=4000;

%--- allocate memory for transition matrix
strans=zeros(4,5); 


    % Load QuB *.dwt file - optional.
    [dwtfilename, dwtpathname]=uigetfile('*.dwt','Select QuB dwt file (optional).');
    if dwtfilename~=0
        dwtfile=strcat(dwtpathname, dwtfilename);
        fid=fopen(dwtfile,'r');
    else
        return;
    end

  
    h=waitbar(0,'Calculating transition matrix.');
    
    % Loop through each trace
    for i=1:ntraces
        seq=[];
        dwells=[];
        line=fgetl(fid);
        
        if ~ischar(line)
            break;
        end
    
        if line(1)=='S'
            line=fgetl(fid);
        end

        while line(1)~='S'
            dwells=[dwells; strread(line)];
            line=fgetl(fid);
            if ~ischar(line)
                break;
            end
        end
%         disp('dwells');disp(dwells);
        seq=dwells(:,1); %--- List of states in one trace
    
        for m=2:numel(seq)
            switch seq(m-1)
                case 0                     %--- ZR state      
                    switch seq(m)
                        case 1             %--- CR state  
                            strans(1,2)=strans(1,2)+1;
                        case 2             %--- GA state 
                            strans(1,3)=strans(1,3)+1;
                        case 3             %--- AC and PB state 
                            strans(1,4)=strans(1,4)+1;
                            if dwells(m,2)>=MAX_THREE_DWELL
                                break
                            end
                    end
                case 1
                    switch seq(m)
                        case 0
                            strans(2,1)=strans(2,1)+1;
                        case 2
                            strans(2,3)=strans(2,3)+1;
                        case 3
                            strans(2,4)=strans(2,4)+1;
                            if dwells(m,2)>=MAX_THREE_DWELL
                                break
                            end
                    end
                case 2
                    switch seq(m)
                        case 0
                            strans(3,1)=strans(3,1)+1;
                        case 1
                            strans(3,2)=strans(3,2)+1;
                        case 3
                            strans(3,4)=strans(3,4)+1;
                            if dwells(m,2)>=MAX_THREE_DWELL
                                break
                            end
                    end
                case 3                              
                    switch seq(m)
                        case 0
                            strans(4,1)=strans(4,1)+1;
                        case 1
                            strans(4,2)=strans(4,2)+1;
                        case 2
                            strans(4,3)=strans(4,3)+1;
                            if dwells(m,2)>=MAX_THREE_DWELL
                                break
                            end
                    end
            end
        
        end
    
        waitbar(i/ntraces,h);
    end
    close(h);
    strans=strans';temp=strans;
    for i=1:5
        temp(6-i,:)=strans(i,:);
    end
    strans=temp;

disp('Transition Matrix (from,to):')
disp(strans(2,:));
disp(strans(3,:));
disp(strans(4,:));
disp(strans(5,:));
