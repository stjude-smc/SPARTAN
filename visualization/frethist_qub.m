function fret_hist=frethist_qub()
%--------------------------------------------------------------------------
%frethist_qub.m takes takes ASCII text data from QuB, and forms 2D
%population histograms with FRET on the y-axis, time on the x-axis, and
%population on the z-axis. At the command line, frethist_qub asks for
%instructions regarding post-synchronization of the FRET traces. There
%are two possible answers:
%           
%           'n': Do not post-synchronize. Outputs space-delimited text file
%                with tag "_frethist".
%           'y': Post-synchronize each FRET traces to the time step at
%                which FRET first occurs within some interval set by the
%                "synch" variables, which are defined below. Output is
%                tagged with "_frethist_ps".
%
%Multiple files may be loaded at one time, outputing one histogram 
%called "FREThist.txt". A list file, also generated in QuB can also be
%entered. This will only place certain regions of each trace into the
%histogram. If no list file is entered, the entirety of the traces are
%used.
%
%NOTE: Unlike frethist, if no list file is entered frethist_qub requires
%that the number of traces, and the length of each trace are entered into
%the variables 'ntraces' and 'len' below.
%
%       -JBM, 3/2006
%--------------------------------------------------------------------------

%   Copyright 2007-2015 Cornell University All Rights Reserved.


fret_hist=[];

%-----Number of traces per file:
ntraces=445;
%-----Number of time points per trace:
len=1501;
%-----Bin size for FRET histogram:
bin=0.025;
%-----FRET value interval between which events are post-synchronized:
%-----synch2 > synch1
synch1=0.125;
synch2=10;
%-----Number of frames before the point of post-synchronization at which
%-----histogram begins:
backset=3;

%---Open the corresonding qub data file
[tracefile tracepath]=uigetfile('*.txt','Choose qub data file:');
if tracefile==0
    disp('No traces file selected.')
    return
end

tracefilename=strcat(tracepath,tracefile);
data=dlmread(tracefilename,' ');

outfile=strrep(tracefilename,'qub.txt','qubhist.txt');
disp('Writing histogram to:'), disp(outfile);


%---Open a QuB list file
[listfile listpath]=uigetfile('*.txt','Choose a QuB list file:');
if listfile==0, 
    disp('No selection list file selected.')
    return;
end
    
listfilename=strcat(listpath,listfile);
list=dlmread(listfilename,'-')+1;
nsegs=size(list,1);

s=size(data);

%-----Axes for histogram:
time_axis=1:len;
fret_axis=-0.2:bin:1.2;

%-----Initialize histogram array, setting the time step in the first row,
%-----and the FRET bins in the first column. This is done for import into
%-----Origin.
fret_hist=zeros(length(fret_axis)+1,length(time_axis)+1);
fret_hist(1,1)=0;
fret_hist(1,2:end)=time_axis;
fret_hist(2:end,1)=fret_axis';

ans=input('Postsynchronize FRET traces (y/n)?','s');
switch (ans)
    case 'n'

        for k=1:nsegs
            seg=data(list(k,1):list(k,2),1);
            for n=1:length(seg)
                fret_hist(2:end,n+1)=fret_hist(2:end,n+1)+...
                    hist(seg(n,1),fret_axis)';
            end
        end
        
    case 'y'
        for k=1:nsegs
            seg=data(list(k,1):list(k,2),1);
            tmp=(seg>=synch1)+(seg<=synch2)-1;
            inds=find(tmp);
            if isempty(inds)
                continue
            end
            pspoint=min(inds);
            start=max([1 pspoint-backset]);
            m=2;
            for n=start:length(seg)
                fret_hist(2:end,m)=fret_hist(2:end,m)+...
                    hist(seg(n,1),fret_axis)';
                m=m+1;
            end
        end
        
end

%-----Write output file:
dlmwrite(outfile,fret_hist,' ');

