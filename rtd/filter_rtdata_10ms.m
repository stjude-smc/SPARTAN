function filter_rtdata_10ms(tracesFile,dwtFile,listFile)
% FILTERRTDATA   
% 
%   Removes traces from datasets not in a given selection list.
%   Filtering criteria used in QuB: DwellCount <= mu + 2sigma
%   Outputs auto.txt, .qub.txt, and .DWT files.
%   Ultimately this script should do the filtering itself.
 

% Get filenames from user
if nargin<1
    [datafile,datapath] = uigetfile({'*.txt'},'Choose a traces.txt file:');
    tracesFile = [datapath datafile];
end
if nargin<2,
    [datafile,datapath] = uigetfile({'*.dwt'},'Choose a QuB DWT file:');
    dwtFile = [datapath datafile];
end
if nargin<3,
    [datafile,datapath] = uigetfile({'*.txt'},'Choose a QuB list file:');
    listFile = [datapath datafile];
end

% if datafile==0
%     disp('No files specified, exiting.');
%     return;
% end


% Load traces file
[d,a,f,ids] = loadTraces( tracesFile );
[nTraces,traceLen] = size(d);

% Load idealization file
[dwells,DT,offsets,model] = loadDWT( dwtFile );

% Load filtered selection list:
% I am assuming here that traces can only be removed; the selection
% ranges are never changed.
fid = fopen( listFile,'r' );
data = textscan(fid,'%d - %d');  %start - end
data = cell2mat( data );
fclose(fid);

startTimes = data(:,1);  %from selection list


% Filter and save DWT data (not used)
keepers = ismember( offsets, startTimes );
dwells = dwells(keepers);

saveDWT( 'flt_qub.dwt', dwells, startTimes, model, DT );



% Find the indexes of traces to keep
selectedTraces   = floor(startTimes/traceLen) +1;

% Filter trace data and save it
d = d(selectedTraces,:);
a = a(selectedTraces,:);
f = f(selectedTraces,:);
ids = ids(selectedTraces);

disp( sprintf('Saving %.0f traces...',size(d,1)) );

saveTraces( 'flt.txt', 'txt', d,a,f,ids );
saveTraces( 'flt_qub.txt', 'qub', f );



