function frethists=statehist(dwtfilename,tracefilename,listfilename, ...
                             lw_limit,up_limit)

%---Builds 2-dimensional histogram of initial and final FRET values for
%---each transition in a group of traces. Data must have been idealized in
%---QuB first, and the dwell times saved in a .dwt file.

%---JBM, 12/06

%---Histogram bin size
BIN=0.015; 
%---Time step (ms)
DT=10;
%---Length of a trace (frames)
TIME=2500;
%--- Integration starts above this lower time limit (ms)
if nargin<4,
    lw_limit=input('lower_limit 10....25000 ms:');
end
if nargin<5,
    up_limit=input('upper_limit 10....25000 ms:');
end
%---Length of dwell (in ms) in state MAX_DWELL_STATE, after which the trace is
%---eliminated. Deactivate this by making MAX_DWELL very large.
MAX_DWELL=20000;
MAX_DWELL_STATE=1;

%---Histogram axes
fret_axis=-0.2:BIN:1.2;

%---Allocate some memory for later
fret=[];

%---Molecule counter
sel_mol_no=0;
mol_no=0;

%---Initialize the histograms
frethists=zeros(numel(fret_axis),5);
frethists(:,1)=fret_axis';

files=0;
% while -1
    list=[];

    %---Open the QuB dwt file from idealization
    if nargin<1
        [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file (2stateModel):');
%     if dwtfile==0
%         break
%     else
        dwtfilename=strcat(dwtpath,dwtfile);
    end
    fid=fopen(dwtfilename,'r');
%     end

    %---Open the corresonding qub data file
    if nargin<2
        [tracefile tracepath]=uigetfile('*.txt','Choose qub data file:');
%     if tracefile==0
%         break
%     else
        tracefilename=strcat(tracepath,tracefile);
    end
    data=dlmread(tracefilename,' ');
    files=files+1;
%     end

    %---Open a QuB list file
    if nargin<3
        [listfile listpath]=uigetfile('*.txt','Choose a QuB list file (required):');
        listfilename=strcat(listpath,listfile);
        
        if listfile==0
            len=length(data);
            list=[1:TIME:len; TIME:TIME:len]';
            nsegs=size(list,1);
        else
            list=dlmread(listfilename,'-')+1;
            nsegs=size(list,1);
        end
    else
        list=dlmread(listfilename,'-')+1;
        nsegs=size(list,1);
    end
    
    
        


    %---Read the dwt file one line at a time
    %---(Looking at the dwt file in a text editor will make clear what's going
    %---on here)
    for i=1:nsegs
        dwells=[];

        line=fgetl(fid);
        if ~ischar(line)
            break;
        end

        if line(1)=='S'
            line=fgetl(fid);
        end

        %---Make a list of dwells from one trace at a time
        while line(1)~='S'
            tmp=strread(line);
            dwells=[dwells; tmp];

            line=fgetl(fid);
            if ~ischar(line)
                break;
            end
        end
    
        if isempty(dwells)
            continue
        end
        times=dwells(:,2)/DT;      %unit of times is frames
        states=dwells(:,1);
        ndwells=numel(times);
        %disp('mol');disp(i);disp(dwells);

        %---Convert the lists of initial and final dwell times into lists of
        %---initial and final FRET values
        seg=data(list(i,1):list(i,2),1);
        
        ulimit=up_limit;
        llimit=lw_limit;
        
        for j=1:ndwells
            if j==1
                ti=1;
                tf=times(j,1);
            else
                ti=tf+1;
                tf=ti+times(j,1)-1;
                
            end
            %disp('(ti,tf');disp(ti);disp(tf);
            
            %if and(states(j)==MAX_DWELL_STATE, times(j)>=MAX_DWELL/DT)  %for truncated traces
            %    for t=1:(tf-ti)
            %        if seg(ti+t-1,1)==seg(ti+t,1)
            %            break
            %        end
            %    end
            %    disp('mol');disp(i);disp(t-1);
            %    hst=hist(seg(ti:ti+t-1,1),fret_axis)';
            %    frethists(:,5)=frethists(:,5)+hst;
            %    break                                         
            %else
            if dwells(j,1)==1
                if and(dwells(j,2)>llimit,dwells(j,2)<ulimit)
                    hst=hist(seg(ti+(llimit/DT):tf,1),fret_axis)';
                    frethists(:,5)=frethists(:,5)+hst;
                elseif and(dwells(j,2)>llimit,dwells(j,2)>=ulimit)
                    hst=hist(seg(ti+(llimit/DT):ti+(ulimit/DT),1),fret_axis)';
                    frethists(:,5)=frethists(:,5)+hst;
                    break
                end
                
                if dwells(j,2)<ulimit
                    ulimit=ulimit-dwells(j,2);
                end
                %disp('up-Limit');disp(ulimit);
                
                if dwells(j,2)<llimit
                    llimit=llimit-dwells(j,2);
                end
                %disp('lw-Limit');disp(llimit);
                
            elseif dwells(j,1)==0
                hst=hist(seg(ti:tf,1),fret_axis)';
                frethists(:,2)=frethists(:,2)+hst;
            end
            
                        
            %switch dwells(j,1)
            %    case 0
            %        frethists(:,2)=frethists(:,2)+hst;
            %    case 1
                    %frethists(:,3)=frethists(:,3)+hst;
            %        frethists(:,5)=frethists(:,5)+hst;
                %case 2
                %    frethists(:,4)=frethists(:,4)+hst;
                %case 3
                %   frethists(:,5)=frethists(:,5)+hst;
            %end
            
            %mean_fret=mean(seg(ti:tf,1))
            %switch mean_fret
             %   case and (mean_fret>-0.125, mean_fret<=0.125)
             %       frethists(:,2)=frethists(:,2)+hst;
             %   case and (mean_fret>0.125, mean_fret<=0.25)
             %       frethists(:,3)=frethists(:,3)+hst;
             %   case and (mean_fret>0.25, mean_fret<=0.45)
             %       frethists(:,4)=frethists(:,4)+hst;
             %   case and (mean_fret>0.45, mean_fret<=0.8)
             %       frethists(:,5)=frethists(:,5)+hst;
            %end
            
        end
        
    end
    fclose(fid);
% end

if files>1
    [outfilename outfilepath]=uiputfile('_sthist.txt','Save tdp file as:');
    outfile=strcat(outfilepath,outfilename);
else
    outfile=strrep(tracefilename,'.txt','_sthist-2.txt');
end

dlmwrite(outfile,frethists,' ');
