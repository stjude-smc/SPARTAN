function shist=statehist(dwtfilename, tracefilename, fret_axis)
% STATEHIST  state occupancy FRET histogram
%
%   S = STATEHIST( DWT, DATA, AXIS )   where
%   S is a collection (in cols) of FRET histograms for each state, as
%   calculated from the DWT file.  The first col specifies the bins.
%   DWT is the filename of the idealization file from QuB.
%   DATA is the auto.txt filename containing raw Fluorescence/FRET data.
%   AXIS specifies FRET histogram bin centers.
%

% Created: Mar 6, 2008    Daniel Terry
% Based on JBM's statehist and tdplot
% TODO: work only on first XXX frames.



% Get filenames from user if not passed
if ~exist('tracefilename','var')
     %---Open the QuB dwt file from idealization
    [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file:');
    if dwtfile==0, return;  end
        
    dwtfilename=strcat(dwtpath,dwtfile);

    %---Open the corresonding qub data file
    [tracefile tracepath]=uigetfile('*.txt','Choose qub data file:');
    if tracefile==0, return;  end
    
    tracefilename=strcat(tracepath,tracefile);
end


%---Histogram axes
if ~exist('fret_axis','var'),
    BIN=0.030;
    fret_axis=-0.15:BIN:0.91;
end


%---Initialize the histogram collection (states in cols)
% shist=zeros( numel(fret_axis), 4+1 );
% shist(:,1)=fret_axis;

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
    %for i=1:nsegs
    while 1,
     
        % Skip segment header line
        tdata = textscan(dwt_fid, 'Segment: %d Dwells: %d Sampling(ms): %d Start(ms): %d ClassCount: %d %*[^\n]');
        if numel(tdata{1}) <= 0,  break;  end  %end of DWT file
        
        nStates = tdata{5};
        %---Initialize the histogram collection (states in cols)
        if ~exist('shist','var')
            shist=zeros( numel(fret_axis), nStates+1 );
            shist(:,1)=fret_axis;
        end
        
        DT = tdata{3};
        if DT==24, DT=25; end
        dwt_start = tdata{4}/DT +1;
        
        % Read dwell states and times
        tdata = textscan(dwt_fid, '%d%d');
        states = tdata{1};
        times  = double(tdata{2}/DT);  %in frames
        
        ndwells=numel(times);
        if ndwells==0, continue; end
        
        dwt_end = dwt_start + sum(times)-1;
        t = t+ndwells-1;  %ntrans=ndwells-1
        
        %total_time = total_time + sum( times(states>0) );
        total_time = total_time + sum( times );
        
        
        %---Sort FRET values into bins for assigned state
        ti = cumsum( [1 ; times(1:end-1)] );
        tf = cumsum( times );
        
        seg_data   = data( dwt_start:dwt_end );
        %centers = (fret_axis(1:end-1)+fret_axis(2:end)) / 2;
        
        for j=1:ndwells,
            points = seg_data( ti(j):tf(j) );
            
            newdata = histc(points,fret_axis);
            newdata = reshape(newdata, numel(fret_axis), 1);
            shist(:,states(j)+2) = shist(:,states(j)+2) + newdata;
        end

    end % for each segment in selection list
    fclose(dwt_fid);
    
% end % for each file


% if files>1
%     [outfilename outfilepath]=uiputfile('_shist.txt','Save shist file as:');
%     outfile=strcat(outfilepath,outfilename);
% else
    outfile=strrep(dwtfilename,'.dwt','_shist.txt');
% end


shist(:,2:end) = shist(:,2:end)/total_time;
dlmwrite(outfile,shist,' ');


% NOTE: the t here is all transitions, including to 0-FRET
% fprintf('t=%d  N=%d  t/n=%f\n', [t nsegs t/nsegs]);

