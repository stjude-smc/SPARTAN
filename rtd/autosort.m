function output = autosort( filename, picks )
% Development of an automated method for the discovery of FRET events
%  above chance background fluctuations.

%--- GLOBAL PARAMTER VALUES

showWB = 1;

constants = cascadeConstants;

%-----FRET value interval between which events are post-synchronized:
%-----synch2 > synch1
% synch1=0.14;
% synch2=10;
%----- Length of Baseline (???)
baseline=100;
%-----Number of time steps within which a Cy3 blink accurs:
blink2=10;
%-----Number of frames before the point of post-synchronization
backset=10;

fretTreshold = 0.14;      %minimum FRET value to consider event
blink = 10;                 %combine events seperated by less than this
scoreThreshold = 0.125;    %accept events that are at least this anticorrelated
areaScoreThreshold = 0.03; %reject short, low FRET events (area under event)

fretOnly = false;  % if true, fluorescence intensity information is unavailable.


%% Load trace data from user
if ~exist('filename','var'),
    [datafile,datapath] = uigetfile({'*.traces'},'Choose a traces file:');
    if datafile==0, return; end

    filename = [datapath filesep datafile];
end

% If the datafile is forQuB, ...
if strfind(filename,'.qub.txt'),
    f = load(filename);
    dataLen = length(f);
    
    fretOnly = true;
    
    % Ask user for number of datapoints per trace, which is unknown
    % in the QuB format...
    traceLen = inputdlg('What is trace length of this data?');
    traceLen = str2double(traceLen);
    
    assert( mod(dataLen,traceLen)==0, 'Incorrect trace length...' );
    nTraces  = dataLen/traceLen;
    
    % Resize data to dimensions of regular auto.txt file (one row per trace).
    f = reshape( f, traceLen,nTraces )';

% If the datafile is an auto.txt file...
else
    if ~exist('picks','var')
        data = loadTraces( filename );
    else
        data = loadTraces( filename, picks );
        showWB = 0;
    end
end


[nTraces,traceLen] = size(data.fret);

% Calculate properties of each trace, which are used for filtering...
if ~fretOnly,
    stats = traceStat( data );
end


%% Seperate events for each trace
% events = zeros(0,2);  %start, stop pairs in rows, one per event

sDonor    = zeros(0,traceLen);
sAcceptor = zeros(0,traceLen);
sFRET     = zeros(0,traceLen);
sIDs      = cell(0,1);
hasEvents = zeros(nTraces,1); %1 if trace has at least one accepted event
nTracesWithEvents = 0;

if showWB
    wb = waitbar(0,'Seperating events...');
end

for m=1:nTraces,

    fret = data.fret(m,:);
        
    %===== Isolate individual events

    % Combine events which are only seperately by a short drop
    % below the FRET threshold
    ind  = fret < fretTreshold;
    indf = 1- rleFilter( ind, blink ); %1 during event

    % Seperate events
    [starts,ends] = rle( indf', 1 );
    starts = max(1,starts-1);
    ends = ends+1;
    nEvents = length(starts);
    
    % By default, accept all events.  Below, some are filtered out...
    accept = starts==starts;


    %--- Determine if an event has anticorrelated edges.
    
    % Calculate a score of the magnitude of anti-correlation in
    % fluorescence (minimumal value of zero).
    if ~fretOnly,
        donor = data.donor(m,:);
        acceptor = data.acceptor(m,:);

        meanT = stats(m).t;
        dd = gradient( donor ) / meanT;
        da = gradient( acceptor ) / meanT;
        et = - (dd.*da) *5;
        et(et<0)=0;
        %clear meanT,dd,da;

        score = et > scoreThreshold; %points above anti-correlation threshold.

        % Accept an event if anti-correlated beyond threshold level.
        passedTreshold = ...
                 ( et(starts)+et(starts+1) ) > scoreThreshold | ...
                 ( et(ends)  +et(ends-1)   ) > scoreThreshold;
        accept( ~passedTreshold ) = false;
    

        % Accept traces if they are terminated by photobleaching
        % accept( ends>lt-3 ) = 1;

        % Special case for 1 frame dwells
        f1 = (starts+1)==(ends-1);
        accept(f1) = score(starts(f1)) & score(ends(f1));

    end
    
    %--- Remove events which repeatedly cross the threshold,
    % indicitive of multiple, poor events.
    % -- not implemented!
    
    %--- Assign a score to how much over the threshold this event lies
    areaScore = zeros( nEvents,1 );

    for i=1:nEvents,
        areaScore(i) = sum(  fret( (starts(i)+1):(ends(i)-1) ) - fretTreshold  );
    end

    accept( areaScore<areaScoreThreshold ) = 0;
    
    
    %--- Accept all long, continuous dwells (FIXME)
    % accept( (ends-1)-(starts+1) >= 10) = 1;
%     accept( areaScore>(0.2*15) ) = 1;


    %--- Remove events falling within a double molecule region of the trace
    if ~fretOnly  %if fluorescence information available...
        safeRegion = stats(m).safeRegion;
        accept( starts<=safeRegion ) = 0;
    end
    
    %---- Remove events that start at time=0
    accept( starts<3 ) = 0;
    
    %---- Remove traces where donor blinks during rise in FRET
    % -- these give no useful information about the process of accomidation
    for i=1:nEvents,
        if accept(i) && any( fret( (starts(i)-3):starts(i) ) == 0 )
            accept(i) = 0;
        end
    end
    
    
    
    %===== Save the data into new traces file data
    starts  = starts(accept);
    ends    = ends(accept);
    nEvents = length(starts);
    
    if nEvents>=1,
        nTracesWithEvents = nTracesWithEvents+1;
        hasEvents(m) = 1;
    end
    
    % Save events into seperate traces
    for i=1:nEvents
        % Leave up <backset> points of background before event
        startp = starts(i)-backset;
        startp = max( 1, startp );

        % Leave up to <baseline> points of background after event
        endp  = ends(i)+baseline;

        if i<nEvents,  %have to trim baseline if it goes into next event
            endp  = min( endp, starts(i+1)-1 );
        end
        endp = min(endp, traceLen);  %don't collect data from next trace

        % Copy data from original trace
        fretSave = fret( startp:endp );
        fretSave( fretSave(1:(backset-1)) >= 0.125 ) = 0.124;
        
        %sFRET(end+1, (startp:endp)-startp+1 ) = fret( startp:endp );
        sFRET(end+1, startp:endp ) = fretSave;

        if ~fretOnly
            sDonor(end+1,:)    = donor;
            sAcceptor(end+1,:) = acceptor;
        end
        
        sIDs{end+1}        = [ ids{m} 'e' num2str(i) ];
    end
    
    if showWB,  waitbar(m/nTraces,wb);  end
    
end %for each trace...

% If no fluorescence data is available, fill it in with zeros
if fretOnly
    sDonor    = zeros( size(sFRET) );
    sAcceptor = zeros( size(sFRET) );
end


if showWB,  close(wb);  end

nEvents = size(sFRET,1);
%disp( sprintf('Found %.0f events from %.0f traces (%.0f total)', ...
%              nEvents, nTracesWithEvents, nTraces) );


%% Save the resulting trace data

% If no ouput arguments specified, save the data to file
if nargout==0
    [datapath] = fileparts(filename);
    if ~isempty(datapath)
        datapath = [datapath filesep];
    end

    % Save each event as a seperate trace
    sData.donor    = sDonor;
    sData.acceptor = sAcceptor;
    sData.fret     = sFRET;
    sData.ids      = sIDs;
    
    saveTraces( [datapath 'ips.traces'], 'traces', sData );
    saveTraces( [datapath 'ips.qub.txt'], 'qub', sFRET );
    
    % Save the set of traces that have at least one event.
    idxTracesWithEvents = find( hasEvents );
    saveTraces( [datapath 'tracesWithEvents.txt'],'txt', data.donor(idxTracesWithEvents,:), ...
        data.acceptor(idxTracesWithEvents,:), data.fret(idxTracesWithEvents,:), ...
        data.ids(idxTracesWithEvents) );
    
% Otherwise, Save FRET data to output
else
    output = sFRET;
end








