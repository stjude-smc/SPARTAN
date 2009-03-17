function tdp=tdplot_rtdata_10ms

%---Builds 2-dimensional histogram of initial and final FRET values for
%---each transition in a group of traces. Data must have been idealized in
%---QuB first, and the dwell times saved in a .dwt file.

%---JBM, 12/06
%---modiefied for rt-Data

%---Histogram bin size
BIN=0.035; 
%---Time step (ms)
DT=10;
%---Length of a trace (frames)
TIME=2500;
%---Length of dwell (in ms) in state MAX_DWELL_STATE, after which the trace is
%---eliminated. Deactivate this by making MAX_DWELL very large.
MAX_DWELL=10000;
MAX_DWELL_STATE=3;

%---Histogram axes
fret_i_axis=-0.2:BIN:1.2;
fret_f_axis=-0.2:BIN:1.2;

%---Allocate some memory for later
h_i=[];
h_f=[];
fret=[];
fret_all=[];
all_rtfret=[];
all_qubfret=[];
time_axis=1:TIME;
rt_fret=zeros(1,TIME);
rt_cy3=zeros(1,TIME);
rt_cy5=zeros(1,TIME);
time_axis=1:TIME;

%---Molecule counter
sel_mol_no=0;
mol_no=0;
%---Initialize the 2d histogram
tdp=zeros(numel(fret_f_axis)+1,numel(fret_i_axis)+1);
tdp(1,2:end)=fret_i_axis;
tdp(2:end,1)=fret_f_axis;

files=0;lst=0;
while 1
    list=[];

    %---Open the QuB dwt file from idealization
    [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file:');
    if dwtfile==0
        break
    else
        dwtfilename=strcat(dwtpath,dwtfile);
        fid=fopen(dwtfilename,'r');
    end

    %---Open the corresonding qub data file
    [tracefile tracepath]=uigetfile('*.txt','Choose qub data file:');
    if tracefile==0
        break
    else
        tracefilename=strcat(tracepath,tracefile);
        data=dlmread(tracefilename,' ');
        files=files+1;
    end

    %---Open a QuB list file
    nodata=[];noempty=[];
    [listfile listpath]=uigetfile('*lst.txt','Choose a QuB list file:');
    if listfile~=0
        listfilename=strcat(listpath,listfile);
        list=dlmread(listfilename,'-')+1;
        list(:,2)=list(:,2)+1;
        nsegs=size(list,1);
        for j=1:nsegs
            nodata=[nodata;((list(j,1)-1)/TIME)+1];
        end
    else
        len=length(data);
        nsegs=(len/TIME);
        list=[1:TIME:len; TIME:TIME:len]';
    end
    
    %---Read the dwt file one line at a time
    %---(Looking at the dwt file in a text editor will make clear what's going
    %---on here)
    empty_trace=0;
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
        seg=data(list(i,1):list(i,2),1);
%         disp('molecule');disp(((list(i,1)-1)/TIME)+1);
%         disp('dwells');disp(dwells);
        
        nonzero=find(dwells(:,2));
        dwells=dwells(size(dwells,1)-size(nonzero,1)+1:end,1:2);
        times=dwells(:,2)/DT; %unit dwells: multiple of 40ms, unit of times is images!
        states=dwells(:,1);
        ndwells=numel(times);
        
        %---
        %---Store Trace as inputfile for FREThist & QUB if not empty
        %---
        rt_fret=zeros(1,TIME);
        qub_fret=zeros(TIME,1);
        if ndwells>1
            rt_fret=seg';
            all_rtfret=[all_rtfret;rt_cy3;rt_cy5;rt_fret];
            all_qubfret=[all_qubfret;seg];
        else 
            empty_trace=empty_trace+1;
        end
        
        %---
        %---Convert the lists of initial and final dwell times into lists of
        %---initial and final FRET values
        %---
        for j=1:ndwells
            if states(j)==MAX_DWELL_STATE & times(j)>=MAX_DWELL
                break
            end
            if j==1 
                ti=1;
                tf=times(j,1);
            else
                ti=tf+1;
                tf=ti+times(j,1)-1;
            end
            
            if states(j)==MAX_DWELL_STATE & times(j)>=3
                fret=[fret mean(seg(ti:ti+3,1))];
            else
                fret=[fret mean(seg(ti:tf,1))];
            end
            fret_all=[fret_all; mean(seg(ti:tf,1))];
        end
     
        ndwells=numel(fret);
        
        %---Add fret values to histogram.
        for k=1:ndwells-1
            h_i=hist(fret(k),fret_i_axis);
            h_f=hist(fret(k+1),fret_f_axis);
            indi=find(h_i>0);
            indf=find(h_f>0);
            tdp(indf+1,indi+1)=tdp(indf+1,indi+1)+1;
        end
        fret=[];
    end
    fclose(fid);
end
ss=sum(sum(tdp(2:end,2:end)));
disp('Total Number of Transitions:'), disp(ss(1));
disp('No. of empty traces');disp(empty_trace);
disp('No. of data traces');disp(size(all_rtfret,1)/3);


%--- Normalize Histogram
ans=input('Constant Normalization Factor (y/n)?','s');
switch (ans)

    %-----Normalization of TDplot:
    case 'n'
        tdp=100*(tdp/ss(1));
        tdp(1,2:end)=0;
        tdp(2:end,1)=0;
        summe=sum(sum(tdp));
        disp('sum of histogram = ');disp(summe);
        disp('Norm.Factor = ');disp(ss(1));
        tdp(1,2:end)=fret_i_axis;
        tdp(2:end,1)=fret_f_axis;
    case 'y'
        norm_factor=input('Normalization Factor:');
        tdp=100*(tdp/norm_factor);
        tdp(1,2:end)=0;
        tdp(2:end,1)=0;
        summe=sum(sum(tdp));
        disp('sum of histogram = ');disp(summe);
        tdp(1,2:end)=fret_i_axis;
        tdp(2:end,1)=fret_f_axis;
end



%--- Save Files
if files>1
    [outfilename outfilepath]=uiputfile('_tdp.txt','Save tdp file as:');
    outfile=strcat(outfilepath,outfilename);
else
    outfile=strrep(tracefilename,'.txt','_tdp.txt');
end
dlmwrite(outfile,tdp,' ');
%dlmwrite('all_fret.txt',fret_all,' ');
%dlmwrite('_flt-ips.txt',all_rtfret,' ');
%dlmwrite('_flt-ips_qub.txt',all_qubfret,' ');
