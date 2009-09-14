function tdp=tdplot(dwtfilename,tracefilename,varargin)

%---Builds 2-dimensional histogram of initial and final FRET values for
%---each transition in a group of traces. Data must have been idealized in
%---QuB first, and the dwell times saved in a .dwt file.

% Performance note: textscan > load > dlmread for speed




% Get filenames from user if not passed
if nargin<2
     %---Open the QuB dwt file from idealization
    [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file:');
    if dwtfile==0, return;  end
        
    dwtfilename=strcat(dwtpath,dwtfile);

    %---Open the corresonding qub data file
    [tracefile tracepath]=uigetfile('*.txt','Choose an auto.txt data file:');
    if tracefile==0, return;  end
    
    tracefilename=strcat(tracepath,tracefile);
end

% Load default values for plotting options.
constants = cascadeConstants;
options.tdp_fret_axis = constants.tdp_fret_axis;
options.normalize = 'total time';

% Modify options if they are specified in the argument list.
if nargin==2,
    assert( isstruct(varargin{1}) );
    options = catstruct( options, varargin{1} );

elseif nargin>3,    
    assert( iscell(varargin) & mod(numel(varargin),2)==0, ...
            'Incorrect format for optional arguments list' );
    vopt = struct(varargin{:});
    options = catstruct( options, vopt );
end




%%

%---Histogram axes
fret_axis = options.tdp_fret_axis;


%---Initialize the 2d histogram
tdp=zeros( numel(fret_axis)+1 );
tdp(1,2:end)=fret_axis;
tdp(2:end,1)=fret_axis;

nTrans = 0;   %total number of transitions
total_time = 0;  %total time in frames

%---Open the QuB dwt file from idealization
[dwt,DT,offsets] = loadDWT(dwtfilename);
nTraces = numel(dwt);

%---Open the corresonding qub data file (slowest step)
[d,a,data] = loadTraces(tracefilename);
data = data';
data = data(:);
clear d; clear a;



%---Count all transitions and add to a transition-density matrix
for i=1:nTraces
    
    % Read dwell states, times, and FRET data
    states = double( dwt{i}(:,1) )-1;
    times  = double( dwt{i}(:,2) );  % dwell times in frames
      
    % Get FRET data corrosponding to this trace.
    dwt_start = offsets(i)+1;
    dwt_end   = dwt_start + sum(times)-1;
    
    % Verify that DWT data bounds make sense.
    assert(dwt_end<=numel(data),'DWT index exceeds data size. Are you sure framerate is correct? (%d)', DT);
    fretData  = data( dwt_start:dwt_end );
    
    % If the option is specified, ignore dark state dwells (see makeplots.m)
    if isfield(options,'hideBlinksInTDPlots') && options.hideBlinksInTDPlots,
        
        % Find dark state dwell regions in the FRET data and remove them.
        idl = dwtToIdl( dwt(i), sum(times), 0 );
        fretData = fretData( idl>1 );
        
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
    
    ndwells=numel(times);
    if ndwells==0, continue; end
    
    % Make sure there are no duplicate states.
    assert( all(diff(states)~=0 ),'Duplicate states!' );
    
    % Save data for calculating transitions per second.
    nTrans = nTrans+ndwells-1;  %ntrans=ndwells-1
    total_time = total_time + sum( times(states>0) );

    
    %---Convert the lists of initial and final dwell times into lists of
    %---initial and final FRET values (mean over dwell)
    ti = cumsum( [1 ; times(1:end-1)] );
    tf = cumsum( times );
    
    fret=zeros(ndwells,1);  % mean FRET value of each dwell (1xN)
    for j=1:ndwells
        fret(j) = mean(  fretData( ti(j):tf(j) )  );
    end

    
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
        tdp(indf,indi) = tdp(indf,indi)+1;
    end

end % for each segment in selection list



% Normalize the plot
tdpData = tdp(2:end,2:end);

if strcmpi(options.normalize,'total time')
    maxVal = total_time*DT/1000;  %in sec -- independant of framerate

    tdp(1,1) = maxVal;  %save the normalization factor (not plotted)
    tdp(2:end,2:end) = tdpData/double(maxVal);
    
elseif strcmpi(options.normalize,'total transitions')
    maxVal = sum( tdpData(:) )/100; %percent in each bin
    
    tdp(1,1) = maxVal; %max value=total # of transitions
    tdp(2:end,2:end) = tdpData/double(maxVal);
    
end

% Save the normalized plot to disk.
outfile=strrep(dwtfilename,'.dwt','_tdp.txt');
dlmwrite(outfile,tdp,' ');


% NOTE: the t here is all transitions, including to 0-FRET
fprintf('t=%d  N=%d  t/n=%f\n', [nTrans nTraces nTrans/nTraces]);

