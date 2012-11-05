function [nTraces,nTrans,tdp]=tdplot_rtdata_v3(dwtfilename,tracefilename,normMethod,norm)

%---Builds 2-dimensional histogram of initial and final FRET values for
%---each transition in a group of traces. Data must have been idealized in
%---QuB first, and the dwell times saved in a .dwt file.

%---JBM, 12/06
%---modified for rt-Data
%---DT: modified to simplify code for automation (9/14/09)


%---Length of dwell (in ms) in state MAX_DWELL_STATE, after which the trace is
%---eliminated. Deactivate this by making MAX_DWELL very large.
MAX_DWELL=10000;
MAX_DWELL_STATE=4;

%---Histogram axes
BIN=0.035;
fret_i_axis=-0.2:BIN:1.2;
fret_f_axis=fret_i_axis;


%---Initialize the 2d histogram
tdp=zeros(numel(fret_f_axis)+1,numel(fret_i_axis)+1);


%---Open the QuB dwt file from idealization
if nargin<1,
    [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file:');
    if dwtfile==0
        return;
    end
    
    dwtfilename=strcat(dwtpath,dwtfile);
end

%---Open the corresonding FRET data file
if nargin<2,
    [tracefile tracepath]=uigetfile('*.txt','Choose FRET data file:');
    if tracefile==0
        return;
    end
    
    tracefilename=strcat(tracepath,tracefile);
end


% Load FRET data.
data = loadTraces(tracefilename);
fretData = data.fret;
[nTraces,traceLen] = size(fretData);
fretData = fretData';
fretData = fretData(:);
% Load dwell-time information
[dwt,sampling,offsets] = loadDWT(dwtfilename);
nsegs = numel(dwt);


%---Read the dwt file one line at a time
for i=1:nsegs
    
    % Load dwell-time information for this trace.
    dwells = dwt{i};
    
    if isempty(dwells), continue; end

    times=dwells(:,2); %unit dwells: multiple of 10ms, unit of times is images!
    states=dwells(:,1); %starts at 1, unlike DWT file.
    ndwells=numel(times);

    
    % Load FRET data for this trace.
    seg = fretData( offsets(i)+(1:traceLen-1) );
 
    %---
    %---Convert the lists of initial and final dwell times into lists of
    %---initial and final FRET values
    %---
    fret=[];
    %disp(i);
    %disp('states');disp(states);
    for j=1:ndwells
        if states(j)==MAX_DWELL_STATE && times(j)>=MAX_DWELL
            break
        end
        if j==1 
            ti=1;
            tf=times(j,1);
        else
            ti=tf+1;
            tf=ti+times(j,1)-1;
        end
        
        if tf <traceLen
            fret=[fret mean(seg(ti:tf,1))];
        else 
            fret=[fret mean(seg(ti:traceLen-1,1))];
        end
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
end

ss=sum(sum(tdp(2:end,2:end)));
nTrans=ss(1);

%--- Normalize Histogram
if nargin>=2,
    answer = num2str(normMethod);
else
    answer=input('Transitions (1), Constant (2), Not_Normalized (3)?','s');
end
switch (answer)

    %-----Normalization of TDplot:
    case '1'
        tdp=100*(tdp/ss(1));
        tdp(1,2:end)=0;
        tdp(2:end,1)=0;
        summe=sum(sum(tdp));
        %disp('sum of histogram = ');disp(summe);
        %disp('Norm.Factor = ');disp(ss(1));
        tdp(1,2:end)=fret_i_axis;
        tdp(2:end,1)=fret_f_axis;
    case '2'
        %norm_factor=input('Normalization Factor:');
        norm_factor=norm;
        tdp=100*(tdp/norm_factor);
        tdp(1,2:end)=0;
        tdp(2:end,1)=0;
        summe=sum(sum(tdp));
        %disp('sum of histogram = ');disp(summe);
        tdp(1,2:end)=fret_i_axis;
        tdp(2:end,1)=fret_f_axis;
    case '3'
        %disp('Norm.Factor = ');disp(ss(1));
        tdp(1,2:end)=fret_i_axis;
        tdp(2:end,1)=fret_f_axis;
end

%--- Save Files
if normMethod==1
    outfile=strrep(tracefilename,'.txt','_tdp1.txt');
elseif normMethod==2
    outfile=strrep(tracefilename,'.txt','_tdp2.txt');
elseif normMethod==3
    outfile=strrep(tracefilename,'.txt','_tdp3.txt'); 
end
dlmwrite(outfile,tdp,' ');



