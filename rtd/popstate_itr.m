function dwthst=popstate_itr(dwtfilename)

%BIN=10; %in ms
BIN=1; %in ms
LEN=25000; %in ms 
TIME=2500;
DT=10;

%--- time to peptide bond formation (multiple of BIN)
popzero=zeros(1,LEN);
popone=zeros(1,LEN);
poptwo=zeros(1,LEN);
popthree=zeros(1,LEN);

%-----Number of frames before the point of post-synchronization at which
%-----histogram begins:
backset=80; % in ms

%-----Initialze molecule counter 
mols=0;
all_events=0;

% Load dwell-time information
[dwt,sampling,offsets] = loadDWT(dwtfilename);
nsegs = numel(dwt);

%---Read the dwt file one line at a time
for i=1:nsegs
    
    % Load dwell-time information for this trace.
    dwells = dwt{i};
    
    if isempty(dwells), continue; end

    times=dwells(:,2)*DT; %unit dwells: multiple of 10ms, unit of times is images!
    states=dwells(:,1)-1; %starts at 1, unlike DWT file.
    ndwells=numel(times);
    
    if states(1)>0
        continue
    end
    n=1;event=0;z=0;sum_dwl=0;
    for j=1:ndwells
        dwl=times(j);
        %--- Post-synchronize to the first idl-event:    
        if ndwells==1
            continue
        elseif n==1
            j1=1;
            j2=backset; %---the firtst dwell is set to backset
            mols=mols+1;
        else
            j1=j2+1;
            j2=j1+dwl-1;
            sum_dwl=sum_dwl+dwl;
            if j==ndwells
                j2=LEN;
            end
        end
                
        switch states(j)
            case 0
                popzero(1,j1:j2)=popzero(1,j1:j2)+1;
            case 1
                popone(1,j1:j2)=popone(1,j1:j2)+1;
            case 2
                poptwo(1,j1:j2)=poptwo(1,j1:j2)+1;
            case 3
                popthree(1,j1:j2)=popthree(1,j1:j2)+1;
        end
            n=n+1;
    end
        %disp('event/trace');disp(event);
        all_events=all_events+event;
end

%------Time
Time=((1:TIME*DT)-backset);
Time=Time'/1000;
Time=Time(10:BIN:LEN);
%-----Plot Data and write output file:
popzero=popzero'/mols; popzero=popzero(10:BIN:LEN);
popone=popone'/mols;popone=popone(10:BIN:LEN);
poptwo=poptwo'/mols;poptwo=poptwo(10:BIN:LEN);
popthree=popthree'/mols;popthree=popthree(10:BIN:LEN);
dwthst=[Time,popzero,popone,poptwo,popthree];

%disp('No. of analyzed molecules');disp(mols);

%--- Save Files
outfile=strrep(dwtfilename,'qub.dwt','dwthst.txt');
dlmwrite(outfile,dwthst,' ');