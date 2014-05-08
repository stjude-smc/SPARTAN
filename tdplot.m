function tdp=tdplot(dwtfilename,traces_input,varargin)
%
%
%
%
% NOTE: normalization for total time does not include time spent in dark
% states (blinking)! Transitions to and from the dark state /are/ included,
% however.

%---Builds 2-dimensional histogram of initial and final FRET values for
%---each transition in a group of traces. Data must have been idealized in
%---QuB first, and the dwell times saved in a .dwt file.

% Performance note: textscan > load > dlmread for speed




% Get filenames from user if not passed
if nargin<2
     %---Open the QuB dwt file from idealization
    [dwtfile dwtpath]=uigetfile('*.dwt','Choose QuB dwt file:');
    if dwtfile==0, return;  end
        
    dwtfilename = strcat(dwtpath,dwtfile);

    %---Open the corresonding traces file
    [tracefile tracepath]=uigetfile('*.traces','Choose an traces file:');
    if tracefile==0, return;  end
    
    traces_input = strcat(tracepath,tracefile);
end


% Load dwell-time information (dwt)
[dwt,DT,offsets] = loadDWT(dwtfilename);


% Load FRET data
if ischar(traces_input),
    d = loadTraces(traces_input);
    fret = d.fret;
    
elseif isstruct(traces_input)
    fret = traces_input.fret;
    
elseif isnumeric(traces_input)
    fret = traces_input;
    
else
    error('tdplot: Invalid traces input');
end

data = fret';
data = data(:);



% Load default values for plotting options.
constants = cascadeConstants;
options = constants.defaultMakeplotsOptions;
options.normalize = 'total time';

% Modify options if they are specified in the argument list.
if nargin==3,
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
fret_axis = options.fret_axis;


%---Initialize the 2d histogram
tdp=zeros( numel(fret_axis)+1 );
tdp(1,2:end)=fret_axis;
tdp(2:end,1)=fret_axis;

nTrans = 0;   %total number of transitions
total_time = 0;  %total time in frames

nTraces = numel(dwt);

% Truncate the idealization to match the contour plots, which only show
% some of the data (e.g., first 50 frames).
truncate = false;
if isfield(options,'truncate_tdplot'),
    truncate = options.truncate_tdplot;
end


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
        indsToSave = find( states>0 );
        states = states(indsToSave);
        times  = times(indsToSave);
        
        if isempty(indsToSave), continue; end
        
        % Combine dwells that are now in the same state by converting into
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

    
    %---Convert the lists of initial and final dwell times into lists of
    %---initial and final FRET values (mean over dwell)
    ti = cumsum( [1 ; times(1:end-1)] );
    tf = cumsum( times );
    
    % Only consider FRET data within the user-specified range.
    if truncate,
        keep = ti<options.contour_length;
        ti = ti(keep);  tf = tf(keep);
        tf(end) = min( tf(end), options.contour_length );  %truncate last dwell if necessary.
        ndwells = numel(ti);
    end
    
    fret = zeros(ndwells,1);  % mean FRET value of each dwell (1xN)
    for j=1:ndwells,
        fret(j) = mean(  fretData( ti(j):tf(j) )  );
    end
    
    % Save data for calculating transitions per second.
    % NOTE: this code originally did not consider blinking times in the
    % calculation of total time!
    nTrans = nTrans + ndwells-1;  %ntrans=ndwells-1
    total_time = total_time + tf(end);

    
    % Place FRET values into contour bins.
    inds = zeros( ndwells,1 );  %fret bin number assignment of each dwell
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

else
    error('Invalid normalization setting');
end

% Save the normalized plot to disk.
if nargout==0,
    outfile=strrep(dwtfilename,'.dwt','_tdp.txt');
    dlmwrite(outfile,tdp,' ');
end

% NOTE: the t here is all transitions, including to 0-FRET
fprintf('t=%d  N=%d  t/n=%f\n', [nTrans nTraces nTrans/nTraces]);

