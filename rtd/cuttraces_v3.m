function [time2Peptide, NoCutEvent]=cuttraces_v3(traceFilename,dwtFilename,FRETac,time2pep,sampling)


%---- USER TUNABLE PARAMETERS.
%---State which will be truncated
cutoffState=3;
%--- time2Peptide: the trace is truncated in high fret after being more
%--- than 120ms in high fret. time2Peptide must be a multiple of 40ms  
time2Peptide=time2pep;
%--- For calc. of the mean Fret value after a trace is truncated.  
highFRET=FRETac;
std_highFRET=0.061;

% Level (standard deviation) of noise to add to AC120-noise data
% after accommodation.
sig=0.061;

%Allocation of Memory
fret_i=[];
fret_f=[];

%---Molecule counter
sel_mol_no=0;

%---Open the corresonding file of traces
if nargin<1,
    [filename,filepath]=uigetfile('*.txt','Choose a traces (.txt) file:');
    if filename==0
        disp('No traces file selected.');
        return;
    else
        traceFilename=strcat(filepath,filename);
    end
end

[cy3,cy5,fret,ids,time_axis] = loadTraces(traceFilename);
[nTraces,traceLen] = size(cy3);
fret_org = fret; %fret unmodified
fret_noise = fret; %fret with noise
fret_noNoise = fret; %fret without noise

%---Open the QuB dwt file from idealization
if nargin<1
    [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file:');
    if dwtfile==0
        disp('No dwt file selected.')
        return;
    else
        dwtFilename=strcat(dwtpath,dwtfile);
    end
end

[dwt,DT] = loadDWT(dwtFilename);
assert( numel(dwt)==nTraces, 'Mismatch between FRET data and idealization' );
assert( DT==sampling, 'Unexptected time resolution' );

%---Read the dwt file one line at a time
%---(Looking at the dwt file in a text editor will make clear what's going
%---on here)
emptyTrace=0; allDwells=0;
cutEvent=0;cutMol=0;
pepSelected   = []; %indexes of molecules that form a peptide bond
noPepSelected = []; %indexes of molecules that do not.


wbh = waitbar(0,'Cutting traces...');

for mol_no=1:nTraces,
    dwells = dwt{mol_no};
    dwells(:,1) = dwells(:,1)-1;
    dwells(:,2) = dwells(:,2)*DT;
    
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
    
    ndwells=size(dwells,1);       %number of dwells in the current trace
    state_dwells=dwells;
    dwelltimes=dwells(:,2)/DT;  %--- dwelltimes = 2dim array (1st col=FRETstate, 2nd col= dwell in that state) here select dwell and set unit to 1. 
    
    %--- count empty traces
    if ndwells<=1,
       emptyTrace=emptyTrace+1;
    end
    %--- count overall no of dwells
    allDwells=allDwells+ndwells;
    
    %---Convert the lists of initial and final dwell times into lists of
    %---initial and final FRET values
    sel_mol_no=sel_mol_no+1;
    if ndwells>1           %---ndwell=1, trace-ampl is not changing
        for i=1:ndwells-1       
            if i==1
                id=1;                                    
                fd=dwelltimes(i,1);                         
                fret_state1=mean(fret(sel_mol_no,id:fd));
                fret_i=[fret_i; fret_state1];

                id=dwelltimes(i,1)+1;                        
                fd=id+dwelltimes(i+1,1)-1;
                fret_state2=mean(fret(sel_mol_no,id:fd));
                fret_f=[fret_f; fret_state2];
            else   
                fret_i=[fret_i; mean(fret(sel_mol_no,id:fd))];

                id=fd+1;
                fd=id+dwelltimes(i+1,1)-1;
                fret_f=[fret_f; mean(fret(sel_mol_no,id:fd))];
            end
        end
    end

    %--- store states, dwelltime, amplitudes and molecule-no. in one array
    state_dwells_ampl = [ state_dwells zeros(ndwells,5) ];
    if ndwells==1
        state_dwells_ampl(1,3)=mean(fret(mol_no,1:dwelltimes(1,1)));
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

    %---
    %--- cut traces after 120ms in AC-state
    %---
    z=1;time2acc=state_dwells_ampl(z,2);
    rt_fret=zeros(1,size(fret,2));
    rt_fret_noise=zeros(1,size(fret,2));

    for z=2:ndwells,
        if size(state_dwells_ampl,1)==1, break; end
        time2acc=time2acc+state_dwells_ampl(z,2);
        if z==ndwells && state_dwells_ampl(z,1)<cutoffState 
            break;
        elseif state_dwells_ampl(z,1)==cutoffState && ...
               state_dwells_ampl(z,2)>time2Peptide,
            cutEvent=cutEvent+1;
            cutMol=mol_no;
            break;
        end
    end

    cutoff_no=time2acc/DT; %time (frame) when accommodation occurred.

    

    %--- Cut off FRET data after accommodation event.
    % The FRET value after accommodation is set to a fixed value (rt_fret)
    % or to the value with noise added (rt_fret_noise).
    last_dwl=state_dwells_ampl(z,2);
    len_mean=floor(0.5*(last_dwl/DT));
    idl_line=mean(fret(sel_mol_no,cutoff_no-len_mean:cutoff_no));
    if or(idl_line<=highFRET-std_highFRET,idl_line>=highFRET+std_highFRET)
        idl_line=highFRET;
    end

    %--- Store rt_traces in a matrix
    rt_fret(1,1:cutoff_no)=fret(sel_mol_no,1:cutoff_no);
    rt_fret_noise(1,1:cutoff_no)=fret(sel_mol_no,1:cutoff_no);
    
    if cutMol>0
        z=1:(traceLen-(cutoff_no-(last_dwl/DT)+(time2Peptide/DT)));
        noise = idl_line+sig*randn( size(z) );
        
        rt_fret_noise(1,cutoff_no-(last_dwl/DT)+(time2Peptide/DT)+z)=noise;
        rt_fret(1,cutoff_no-(last_dwl/DT)+(time2Peptide/DT)+z)=idl_line;
    end

    %--- Store rt_traces truncated after 120ms in high FRET  
    if cutMol>0, %
        pepSelected(end+1) = sel_mol_no;
        fret(sel_mol_no,:) = rt_fret;
        
    elseif cutMol==0 && ndwells>1,
        noPepSelected(end+1) = sel_mol_no;
    end

    fret_noNoise(mol_no,:) = rt_fret;
    fret_noise(mol_no,:) = rt_fret_noise;
    

    waitbar(mol_no/nTraces,wbh);
    cutMol=0;
    
end %for each trace.


waitbar(0.1,wbh,'Saving Traces...');

if emptyTrace>=1
    disp('no. of empty Traces');disp(emptyTrace);
end
NoCutEvent=100*(cutEvent/mol_no);
disp(['no. Cutoff-Events in %: ', num2str(NoCutEvent)]);
%--- Save Data

%---
%--- all Molecules
%---
%--- all_rtfret: all Molecules, all traces which are being longer 
%--- than 120ms in high FRET are truncated and subsituted by a 
%--- mean mean FRET-value
z = zeros( size(cy3) );
 saveTraces( 'allMol_ac120.txt','txt', z,z, fret_noNoise, ids,time_axis );
 saveTraces( 'allMol_ac120.qub.txt','qub', fret_noNoise );
 
waitbar(0.2,wbh,'Saving Traces...');

 %saveTraces( 'allMol_ac120-noise.txt','txt', z,z, fret_noise, ids,time_axis );
 %saveTraces( 'allMol_ac120-noise.qub.txt','qub', fret_noise );
 

%---
%--- only Molecules forming a peptide bond
%---
waitbar(0.3,wbh,'Saving Traces...');
 saveTraces('PEP120org.txt','txt',cy3(pepSelected,:),cy5(pepSelected,:), ...
            fret_org(pepSelected,:), ids(pepSelected), time_axis );
waitbar(0.4,wbh,'Saving Traces...');
 saveTraces('PEP120.txt','txt',cy3(pepSelected,:),cy5(pepSelected,:), ...
            fret(pepSelected,:), ids(pepSelected), time_axis );
waitbar(0.5,wbh,'Saving Traces...');
 saveTraces( 'PEP120.qub.txt','qub', fret(pepSelected,:) );
 

%---
%--- only Molecules forming no peptide bond
%---
waitbar(0.6,wbh,'Saving Traces...');
 saveTraces('noPep120.txt','txt',cy3(noPepSelected,:),cy5(noPepSelected,:), ...
            fret(noPepSelected,:),ids(noPepSelected), time_axis );
waitbar(0.7,wbh,'Saving Traces...');
 saveTraces( 'noPep120.qub.txt','qub', fret(noPepSelected,:) );
 
waitbar(0.95,wbh);
close(wbh);









