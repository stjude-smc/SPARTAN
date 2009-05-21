function ctr=cuttraces_2_10ms(traceFilename,dwtFilename)


%---Histogram bin size
BIN=0.035; 
%---Time step (ms)
DT=10;
%---Length of a trace (frames)
LenTrace=2500;
%LenTrace=6004;
%---State which will be truncated
cutoffState=3;
cutoffState_GC=2;
%--- time2Peptide: the trace is truncated in high fret after being more
%--- than 120ms in high fret. time2Peptide must be a multiple of 40ms  
time2Peptide=120; 
time2GTP=320;
%--- time2lngPeptide: is the time a molecule needs to form a peptide bond. 
%--- similar to time2lngPeptide=time2acc>320
%--- time2lngPeptide must be a multiple of 40ms 
time2lngPeptide=320;
dis_time=120;
%--- For calc. of the mean Fret value after a trace is truncated.  
midFRET=0.30;
std_midFRET=0.061;
highFRET=0.55;
std_highFRET=0.061;


%---Allocate some memory for later
h_i=[];
h_f=[];
fret_i=[];
fret_f=[];
kinetic=[];
all_rtfret=[];
all_rtfret_noise=[];
all_rtfret_org=[];
all_qubfret=[];
all_qubfret_noise=[];
all_pepfret_org=[];
all_pepfret_cut=[];
all_noPep_qub=[];
all_idlPEP_qub=[];
all_lngPEP_qub=[];
all_PEPcy3cy5fret=[];
all_lngPEPcy3cy5fret=[];
all_idlPEPcy3cy5fret=[];
all_noPEPcy3cy5fret=[];
all_lngtime2acc=[];
all_cutoff=[];
all_lifetime=[];

allGTP_qub=[];
allGTPcy3cy5fret=[];
all_shtGAcy3cy5fret=[];
all_lngGAcy3cy5fret=[];

all_shtGA_qub=[];
all_lngGA_qub=[];
%---Molecule counter
sel_mol_no=0;
mol_no=0;
rej_no=0;
sel_no=0;

%---Open the corresonding file of traces
files=0;
data=[];

if nargin<1,
    [filename,filepath]=uigetfile('*.txt','Choose a traces (.txt) file:');
    if filename==0
        disp('No traces file selected.');
        return;
    else
        traceFilename=strcat(filepath,filename);
    end
end

fid=fopen(traceFilename,'r');
time=strread(fgetl(fid),'%f')';
% d=dir(traceFilename);
sig=textscan(fid,'%s',1);
sig=sig{:};
if isnan(str2double(sig))
    fseek(fid,-ftell(fid),0);
    time=strread(fgetl(fid),'%f')';
    wb=waitbar(0,'Loading traces...');
    ptr=ftell(fid);
    while ptr~=d.bytes
        line=fgetl(fid);
        id=strread(line,'%s',1);
        data_start=numel(cell2mat(id))+1;
        data=[data; strread(line(data_start:end),'%f')'];
        ptr=ftell(fid);
        waitbar(ftell(fid)/d.bytes,wb);
    end
    fclose(fid);
    close(wb);
else
    fclose(fid);
    new_data=dlmread(file,' ');
    data=[data; new_data(2:end,:)];
end

%---We'll only need the FRET traces
cy3=data(1:3:end,:);
cy5=data(2:3:end,:);
fret=data(3:3:end,:);
cy3zero=zeros(1,LenTrace);
cy5zero=zeros(1,LenTrace);
time_axis=1:LenTrace;

frettrace_no=4:3:size(data,1);
fret_subset=data;


%---Open the QuB dwt file from idealization
if nargin<1
    [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file:');
    [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file:');
    if dwtfile==0
        disp('No dwt file selected.')
        return;
    else
        dwtFilename=strcat(dwtpath,dwtfile);
    end
end

fid=fopen(dwtFilename,'r');


%---Read the dwt file one line at a time
%---(Looking at the dwt file in a text editor will make clear what's going
%---on here)
trace_nos=1:3000;
emptyTrace=0; allDwells=0;
cutEvent=0;cutMol=0;

ans=input('Truncate rt-traces (y/n)?','s');
switch (ans)
    case 'y'
        while 1
            dwells=[];
            line=fgetl(fid);
            if ~ischar(line), break; end
            if line(1)=='S'
                continue;
            else
            %---Make a list of dwells from one trace at a time
                while line(1)~='S'
                    dwells=[dwells; strread(line)];
                    line=fgetl(fid);
                    if ~ischar(line), break; end
                end
            end
            mol_no=mol_no+1;
            disp('molecule');disp(mol_no);
            %disp('dwells');disp(dwells);
            %---- make ini. states with negative meanFRET value to state 0 (delete idl mixing errors) 
            for j=1:size(dwells,1)
                if j==1
                    ti=1;
                    tf=dwells(j,2)/DT;
                else
                    ti=tf+1;
                    tf=ti+(dwells(j,2)/DT)-1;
                end
                if mean(fret(mol_no,ti:tf))<0
                    dwells(j,1)=0;
                end
            end
            %---- combine initial 0 Fret states
            for j=2:size(dwells,1)
                if dwells(j,1)>=1
                    break
                elseif dwells(j,1)==0
                    dwells(j,2)=dwells(j-1,2)+dwells(j,2);
                    dwells(j-1,2)=0;
                end
            end
            nonzero=find(dwells(:,2));
            dwells=dwells(size(dwells,1)-size(nonzero,1)+1:end,1:2);
    
            dwells(:,2)/DT; %--- dwells = 2dim array (1st col=FRETstate, 2nd col= dwell in that state) here select dwell and set unit to 1. 
            state_dwells=dwells;
            dwells=dwells(:,2)/DT;  %--- dwells = 2dim array (1st col=FRETstate, 2nd col= dwell in that state) here select dwell and set unit to 1. 
            ndwells=size(dwells,1); %--- ndwells: number which gives the number of dwells in the current trace
            %---Convert the lists of initial and final dwell times into lists of
            %---initial and final FRET values
            if sum(trace_nos==mol_no)
                sel_mol_no=sel_mol_no+1;
                if ndwells>1                                         %---ndwell=1, trace-ampl is not changing
                    for i=1:ndwells-1       
                        if i==1
                            id=1;                                    
                            fd=dwells(i,1);                         
                            fret_state1=mean(fret(sel_mol_no,id:fd));
                            fret_i=[fret_i; fret_state1];

                            id=dwells(i,1)+1;                        
                            fd=id+dwells(i+1,1)-1;
                            fret_state2=mean(fret(sel_mol_no,id:fd));
                            fret_f=[fret_f; fret_state2];
                        else   
                            fret_i=[fret_i; mean(fret(sel_mol_no,id:fd))];

                            id=fd+1;
                            fd=id+dwells(i+1,1)-1;
                            fret_f=[fret_f; mean(fret(sel_mol_no,id:fd))];
                        end
                    end
                end
            end
    
            %--- store states, dwelltime, amplitudes and molecule-no. in one array
            state_dwells_ampl=zeros(size(dwells,1),5);
            state_dwells_ampl=[state_dwells,state_dwells_ampl];
            if ndwells==1
                state_dwells_ampl(1,3)=mean(fret(mol_no,1:dwells(1,1)));
                state_dwells_ampl(1,4)=sel_mol_no;
            elseif ndwells>1
                for i=1:ndwells
                    if and((size(fret_i,1)+1)==ndwells,i<ndwells)
                        state_dwells_ampl(i,3)=fret_i(i);
                        state_dwells_ampl(i,4)=sel_mol_no;
                    elseif and((size(fret_i,1)+1)==ndwells,i==ndwells)
                        state_dwells_ampl(i,3)=fret_f(size(fret_f,1));
                        state_dwells_ampl(i,4)=sel_mol_no;
                    elseif and((size(fret_i,1)+1)>ndwells, i<ndwells)
                        state_dwells_ampl(i,3)=fret_i(size(fret_i,1)-(ndwells-1)+i);
                        state_dwells_ampl(i,4)=sel_mol_no;
                    elseif and((size(fret_i,1)+1)>ndwells,i==ndwells)
                        state_dwells_ampl(i,3)=fret_f(size(fret_f,1));
                        state_dwells_ampl(i,4)=sel_mol_no;
                    end
                end
            end
            %--- count empty traces
            if ndwells==1;
                emptyTrace=emptyTrace+1;
            end
            %--- count overall no of dwells
            allDwells=allDwells+ndwells;
            
            %---
            %--- cut traces after 120ms in AC-state
            %---
            z=1;time2acc=state_dwells_ampl(z,2);
            rt_fret_org=zeros(1,size(fret,2));
            rt_fret=zeros(1,size(fret,2));
            rt_fret_noise=zeros(1,size(fret,2));
            qub_fret=zeros(size(fret,2),1);
            qub_fret_noise=zeros(size(fret,2),1);
            
            while 1 
                if size(state_dwells_ampl,1)==1,break;
                else z=z+1;end
                time2acc=time2acc+state_dwells_ampl(z,2);
                if and(z==ndwells,state_dwells_ampl(z,1)==0) 
                    break
                elseif and(state_dwells_ampl(z,1)==cutoffState,state_dwells_ampl(z,2)>time2Peptide)
                    cutEvent=cutEvent+1;
                    cutMol=mol_no;
                    break
                end
            end
            
            cutoff_no=time2acc/DT;
            
            %--- time to accomodation
            if z==1                                             %---only 1 dwell
                time2acc=0;
            elseif z==ndwells                                   %---no peptide bond formed
                time2acc=0;
            else                                                %---peptide bond formed
                time2acc=time2acc-state_dwells_ampl(1,2)-state_dwells_ampl(z,2)+time2Peptide;
            end
            
            %---
            %--- Calc. idl-Line and store truncated rt_traces
            %---
            
            %--- Calc. idl-Line
            last_dwl=state_dwells_ampl(z,2);
            len_mean=floor(0.5*(last_dwl/DT));
            idl_line=mean(fret(sel_mol_no,cutoff_no-len_mean:cutoff_no));
            if or(idl_line<=highFRET-std_highFRET,idl_line>=highFRET+std_highFRET)
                idl_line=highFRET;
            end
            
            %--- Store rt_traces in a matrix
            z=0;sig=0.061;
            rt_fret(1,1:cutoff_no)=fret(sel_mol_no,1:cutoff_no);
            rt_fret_noise(1,1:cutoff_no)=fret(sel_mol_no,1:cutoff_no);
            rt_fret_org(1,1:cutoff_no)=fret(sel_mol_no,1:cutoff_no);
            
            if cutMol>0
                for z=1:(LenTrace-(cutoff_no-(last_dwl/DT)+(time2Peptide/DT)))
                    rt_fret_noise(1,cutoff_no-(last_dwl/DT)+(time2Peptide/DT)+z)=idl_line+sig*randn;
                    rt_fret(1,cutoff_no-(last_dwl/DT)+(time2Peptide/DT)+z)=idl_line;
                end
            end
    
            %--- Store rt_traces in a QuB-vector
            z=0;
            qub_fret(1:cutoff_no,1)=fret(sel_mol_no,1:cutoff_no);
            qub_fret_noise(1:cutoff_no,1)=fret(sel_mol_no,1:cutoff_no);
            if cutMol>0
                for z=1:(LenTrace-(cutoff_no-(last_dwl/DT)+(time2Peptide/DT)))
                    qub_fret_noise(cutoff_no-(last_dwl/DT)+(time2Peptide/DT)+z,1)=idl_line+sig*randn;
                    qub_fret(cutoff_no-(last_dwl/DT)+(time2Peptide/DT)+z,1)=idl_line;
                end
            end
            
            %--- Store rt_traces truncated after 120ms in high FRET 
            lngtime2acc=zeros(1,2);
            lngPEPcy3cy5fret=zeros(1,4);
            shtGAcy3cy5fret=zeros(1,4);
            lngGAcy3cy5fret=zeros(1,4);
           
    
            if cutMol>0
                all_cutoff=[all_cutoff;cutMol];
                all_pepfret_cut=[all_pepfret_cut;qub_fret];
                all_pepfret_org=[all_pepfret_org;fret(sel_mol_no,1:LenTrace)];
                all_PEPcy3cy5fret=[all_PEPcy3cy5fret;...
                    cy3(sel_mol_no,:);...
                    cy5(sel_mol_no,:);...
                    rt_fret_noise(1,:)];
                if time2acc>time2lngPeptide
                    all_lngPEP_qub=[all_lngPEP_qub;qub_fret];
                    lngPEPcy3cy5fret=[;...
                    cy3(sel_mol_no,:);...
                    cy5(sel_mol_no,:);...
                    fret(sel_mol_no,:)];
                    all_lngPEPcy3cy5fret=[all_lngPEPcy3cy5fret;lngPEPcy3cy5fret];
                    lngtime2acc=[sel_mol_no time2acc];
                    all_lngtime2acc=[all_lngtime2acc;lngtime2acc];
                elseif time2acc<=time2lngPeptide
                    all_idlPEP_qub=[all_idlPEP_qub;qub_fret];
                    idlPEPcy3cy5fret=[;...
                    cy3(sel_mol_no,:);...
                    cy5(sel_mol_no,:);...
                    fret(sel_mol_no,:)];
                    all_idlPEPcy3cy5fret=[all_idlPEPcy3cy5fret;idlPEPcy3cy5fret];
                end
            elseif and(cutMol==0,ndwells>1)
                all_noPep_qub=[all_noPep_qub;qub_fret];
                all_noPEPcy3cy5fret=[all_noPEPcy3cy5fret;...
                    cy3(sel_mol_no,:);...
                    cy5(sel_mol_no,:);...
                    fret(sel_mol_no,:)];
            end
            
            dwl_state23=0; %sum of dwelltime in state 2 and 3
            for i=1:ndwells
                if or(state_dwells_ampl(i,1)==2,state_dwells_ampl(i,1)==3)
                    dwl_state23=dwl_state23+state_dwells_ampl(i,2);
                end
            end 
            
   
            if dwl_state23<=dis_time
                all_shtGA_qub=[all_shtGA_qub;qub_fret];
                shtGAcy3cy5fret=[;...
                cy3(sel_mol_no,:);...
                cy5(sel_mol_no,:);...
                fret(sel_mol_no,:)];
                all_shtGAcy3cy5fret=[all_shtGAcy3cy5fret;shtGAcy3cy5fret];
            elseif dwl_state23>dis_time
                all_lngGA_qub=[all_lngGA_qub;qub_fret];
                lngGAcy3cy5fret=[;...
                cy3(sel_mol_no,:);...
                cy5(sel_mol_no,:);...
                fret(sel_mol_no,:)];
                all_lngGAcy3cy5fret=[all_lngGAcy3cy5fret;lngGAcy3cy5fret];
            end
    
            cutMol=0;
            
            %---- Store Data 
            kinetic=[kinetic;state_dwells_ampl];
            all_qubfret_noise=[all_qubfret_noise;qub_fret_noise];
            all_qubfret=[all_qubfret;qub_fret];
            all_rtfret=[all_rtfret;cy3zero;cy5zero;rt_fret];
            all_rtfret_noise=[all_rtfret_noise;cy3zero;cy5zero;rt_fret_noise];
            all_rtfret_org=[all_rtfret_org;cy3zero;cy5zero;rt_fret_org];
        end
        
        all_rtfret=[time_axis;all_rtfret];
        all_rtfret_noise=[time_axis;all_rtfret_noise];
        all_rtfret_org=[time_axis;all_rtfret_org];
        all_PEPcy3cy5fret=[time_axis;all_PEPcy3cy5fret];
        all_noPEPcy3cy5fret=[time_axis;all_noPEPcy3cy5fret];
        all_shtGAcy3cy5fret=[time_axis;all_shtGAcy3cy5fret];
        all_lngGAcy3cy5fret=[time_axis;all_lngGAcy3cy5fret];
        

        disp('size time_axis');disp(size(time_axis));
        disp('size lng PEP');disp(size(all_lngPEPcy3cy5fret));
        all_lngPEPcy3cy5fret=[time_axis;all_lngPEPcy3cy5fret];
        all_idlPEPcy3cy5fret=[time_axis;all_idlPEPcy3cy5fret];


        disp('size of FRET matrix');disp(size(all_rtfret));
        disp('no. of empty Traces');disp(emptyTrace);
        disp('no. of allDwells');disp(allDwells);
        disp('no. Cutoff-Events');disp(cutEvent);
        fclose(fid);

        %---
        %--- Save Data
        %---
        %dlmwrite('cutoff_no.txt',all_cutoff,' ');
        %dlmwrite('lngTime2Acc.txt',all_lngtime2acc,' ');
        %dlmwrite('kinetic_all.txt',kinetic,' ');
        %---
        %--- all Molecules
        %---
        %--- all_rtfret_org: all Molecules truncated after being longer 
        %--- than 120ms in high FRET
        %--- all_rtfret: all Molecules, all traces which are being longer 
        %--- than 120ms in high FRET are truncated and subsituted by a 
        %--- mean mean FRET-value
        %dlmwrite('_ac120_org.txt',all_rtfret_org,' ');
         dlmwrite('allMol_ac120.txt',all_rtfret,' ');
         dlmwrite('allMol_ac120-noise.txt',all_rtfret_noise,' ');
         dlmwrite('allMol_ac120_qub.txt',all_qubfret,' ');
         dlmwrite('allMol_ac120-noise_qub.txt',all_qubfret_noise,' ');
        %---
        %--- only Molecules forming a peptide bond
        %---
        %dlmwrite('PEP_ac120_qub.txt',all_pepfret_org,' ');
         dlmwrite('PEP120_qub.txt',all_pepfret_cut,' ');
         dlmwrite('PEP120.txt',all_PEPcy3cy5fret,' ');
        %---
        %--- only Molecules going in slow
        %---
        %dlmwrite('lngPEP320_qub.txt',all_lngPEP_qub,' ');
        %dlmwrite('lngPEP320.txt',all_lngPEPcy3cy5fret,' ');
        %---
        %--- only Molecules going in fast
        %---
        %dlmwrite('idlPep320_qub.txt',all_idlPEP_qub,' ');
        %dlmwrite('idlPep320.txt',all_idlPEPcy3cy5fret,' ');
        %---
        %--- only Molecules forming no peptide bond
        %---
         dlmwrite('noPep120_qub.txt',all_noPep_qub,' ');
         dlmwrite('noPep120.txt',all_noPEPcy3cy5fret,' ');
        %---
        %--- only Molecules dissociating from the GA state
        %---
        %dlmwrite('shrtGA120_qub.txt',all_shtGA_qub,' ');
        %dlmwrite('shrtGA120.txt',all_shtGAcy3cy5fret,' ');
        %disp('shtGA');disp(size(all_shtGAcy3cy5fret));
        %---
        %--- only Molecules with lng lifetime in the GA state
        %---
        %dlmwrite('lngGA120_qub.txt',all_lngGA_qub,' ');
        %dlmwrite('lngGA120.txt',all_lngGAcy3cy5fret,' ');
        %disp('lngGA');disp(size(all_lngGAcy3cy5fret));
        
    case 'n'
        disp('No traces have been modified.');
end
