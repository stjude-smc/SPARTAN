function tdp=tdplot(dwtfilename,tracefilename,options)

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
    options = cascadeConstants();
end


%---Histogram axes
fret_axis = options.tdp_fret_axis;


%---Initialize the 2d histogram
tdp=zeros( numel(fret_axis)+1 );
tdp(1,2:end)=fret_axis;
tdp(2:end,1)=fret_axis;

nTrans = 0;   %total number of transitions
total_time = 0;  %total time in frames

%---Open the QuB dwt file from idealization
[dwt,sampling,offsets] = loadDWT(dwtfilename);
nTraces = numel(dwt);

%---Open the corresonding qub data file (slowest step)
[d,a,data] = loadTraces(tracefilename);
data = data';
data = data(:);
clear d; clear a;


%---Count all transitions and add to a transition-density matrix
for i=1:nTraces

    % Read dwell states, times, and FRET data
    states = double( dwt{i}(:,1) );
    times  = double( dwt{i}(:,2) );  % dwell times in frames
    
    trace = data( offsets(i)+(1:sum(times)) );

    % If the option is specified, ignore dark state dwells (see makeplots.m)
    if isfield(options,'hideBlinksInTDPlots') && options.hideBlinksInTDPlots,
        
        % Find dark state dwell regions in the FRET data and remove them.
        idl = dwtToIdl( dwt(i), sum(times), 0 );
        trace = trace( idl>1 );
        
        % Remove dwells in lowest FRET state (assuming it is the dark state)
        indsToSave = states>1;
        states = states(indsToSave);
        times  = times(indsToSave);
        
        % Combine dwells that are now in the same state by converting into%
        % an idealization and then back to a dwell-time sequence.
        idl = dwtToIdl( {[states times]}, sum(times), 0 );
        dwells = RLEncode(idl);
        states = dwells(:,1);
        times  = dwells(:,2);
    end
    
    assert( all(diff(states)~=0 ) );
    
    nDwells = numel(times);
    if nDwells<1, continue; end

    nTrans = nTrans+nDwells-1;
    total_time = total_time + sum(times);


    %---Convert the lists of initial and final dwell times into lists of
    %---initial and final FRET values (mean over dwell)
    if offsets(i)+sum(times) > numel(data),
        disp(tracefilename);
        error('DWT index exceeds data size. Are you sure the framerate is correct? (%d)', sampling);
    end
    
    
    ti = cumsum( [1 ; times(1:end-1)] );
    tf = cumsum( times );

    fret=zeros(1,nDwells);  % mean FRET value of each dwell (1xN)
    for j=1:nDwells
        fret(j) = mean(  trace( ti(j):tf(j) )  );
    end
    
    %assert(  all( abs(diff(fret))>0.1 )  );
    

    % Place FRET values into contour bins
    inds = zeros( nDwells,1 );  %bin number of each dwell
    centers = (fret_axis(1:end-1)+fret_axis(2:end)) / 2;

    inds( fret<=centers(1) ) = 1;
    for j=1:numel(centers),
        inds( fret>centers(j) ) = j+1;
    end

    % Add transitions to TD plot
    for j=1:nDwells-1            
        indi = inds(j)+1;
        indf = inds(j+1)+1;
        tdp(indf,indi) = tdp(indf,indi)+1;
    end

end % for each segment in selection list


outfile=strrep(dwtfilename,'.dwt','_tdp.txt');

% Normalize the plot and write to disk
max_val = total_time*sampling/1000;  %in sec -- independant of framerate

tdp(1,1) = max_val;  %save the normalization factor (not plotted)
tdp(2:end,2:end) = tdp(2:end,2:end)/double(max_val);
dlmwrite(outfile,tdp,' ');



% NOTE: the t here is all transitions, including to 0-FRET
fprintf('t=%d  N=%d  t/n=%f\n', [nTrans nTraces nTrans/nTraces]);

