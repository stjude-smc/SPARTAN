function tdp=tdplot(dwtfilename,tracefilename,constants)

%---Builds 2-dimensional histogram of initial and final FRET values for
%---each transition in a group of traces. Data must have been idealized in
%---QuB first, and the dwell times saved in a .dwt file.

% Performance note: textscan > load > dlmread for speed

%---JBM, 12/06

% Get filenames from user if not passed
if ~exist('tracefilename','var')
     %---Open the QuB dwt file from idealization
    [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file:');
    if dwtfile==0, return;  end
        
    dwtfilename=strcat(dwtpath,dwtfile);

    %---Open the corresonding qub data file
    [tracefile tracepath]=uigetfile('*.txt','Choose an auto.txt data file:');
    if tracefile==0, return;  end
    
    tracefilename=strcat(tracepath,tracefile);
end

% Load constants used
if nargin<3,
    constants = cascadeConstants();
end


%---Histogram axes
fret_axis = constants.tdp_fret_axis;


%---Initialize the 2d histogram
tdp=zeros( numel(fret_axis)+1 );
tdp(1,2:end)=fret_axis;
tdp(2:end,1)=fret_axis;

t = 0;   %total number of transitions
total_time = 0;  %total time in frames
% files=0;
% while 1
%     list=[];

    %---Open the QuB dwt file from idealization
    dwt_fid=fopen(dwtfilename,'r');

    %---Open the corresonding qub data file (slowest step)
    [d,a,data] = loadTraces(tracefilename);
    data = data';
    data = data(:);
    clear d; clear a;
    
    
%     files=files+1;

    %---Read the dwt file one line at a time
    %---(Looking at the dwt file in a text editor will make clear what's going
    %---on here)
    nsegs = 0;
    %for i=1:nsegs
    while 1,
     
        % Skip segment header line
        tdata = textscan(dwt_fid, 'Segment: %d Dwells: %d Sampling(ms): %d Start(ms): %d %*[^\n]');
        if numel(tdata{1}) <= 0,  break;  end  %end of DWT file
        
        DT = tdata{3};
        if DT==24, DT=25; end
        dwt_start = tdata{4}/DT +1;
        
        % Read dwell states and times
        tdata = textscan(dwt_fid, '%d%d');
        states = tdata{1};
        times  = double(tdata{2}/DT);  %in frames
        
        ndwells=numel(times);
        if ndwells==0, continue; end
        nsegs = nsegs+1;
        
        dwt_end = dwt_start + sum(times)-1;
        t = t+ndwells-1;  %ntrans=ndwells-1
        total_time = total_time + sum( times(states>0) );
        
        
        %---Convert the lists of initial and final dwell times into lists of
        %---initial and final FRET values (mean over dwell)
        if dwt_end>numel(data),
            disp(tracefilename);
            error('DWT index exceeds data size. Are you sure framerate is correct? (%d)', DT);
        end
        seg_data = data( dwt_start:dwt_end );
        fret=zeros(1,ndwells);  % mean FRET value of each dwell (1xN)
        
        
        ti = cumsum( [1 ; times(1:end-1)] );
        tf = cumsum( times );
        
        for j=1:ndwells
%             if tf(j) > size(seg_data,1)  continue;  end
            fret(j) = mean(  seg_data( ti(j):tf(j) )  );
        end
        
        
    
        %---Redefine ndwells to reflect only those dwells before the
        %---MAX_DWELL in MAX_DWELL_STATE
        %ndwells=numel(fret);
        
        %---Add fret values to histogram.
%         for k=1:ndwells-1
%             h_i=hist(fret(k),fret_i_axis);
%             h_f=hist(fret(k+1),fret_f_axis);
%             indi=find(h_i>0);
%             indf=find(h_f>0);
%             tdp(indf+1,indi+1)=tdp(indf+1,indi+1)+1;
%         end
        
        % Place FRET values into contour bins
        inds = zeros( ndwells,1 );  %bin number of each dwell
        centers = (fret_axis(1:end-1)+fret_axis(2:end)) / 2;
        
        inds( fret<=centers(1) ) = 1;
        for n=1:numel(centers),
            inds( fret>centers(n) ) = n+1;
        end
        
        % Add transitions to TD plot
        for k=1:ndwells-1            
            indi = inds(k)+1;
            indf = inds(k+1)+1;
%             tdp(indf,indi)=tdp(indf,indi)+(1/(ndwells-1));
            tdp(indf,indi) = tdp(indf,indi)+1;
        end
        
    end % for each segment in selection list
    fclose(dwt_fid);
    
% end % for each file


% if files>1
%     [outfilename outfilepath]=uiputfile('_tdp.txt','Save tdp file as:');
%     outfile=strcat(outfilepath,outfilename);
% else
    outfile=strrep(dwtfilename,'.dwt','_tdp.txt');
% end

% Normalize the plot and write to disk
% max_val = max(max( tdp(2:end,2:end) ));
max_val = total_time*DT/1000;  %in sec -- independant of framerate

tdp(1,1) = max_val;  %save the normalization factor (not plotted)
tdp(2:end,2:end) = tdp(2:end,2:end)/double(max_val);
dlmwrite(outfile,tdp,' ');



% NOTE: the t here is all transitions, including to 0-FRET
fprintf('t=%d  N=%d  t/n=%f\n', [t nsegs t/nsegs]);

