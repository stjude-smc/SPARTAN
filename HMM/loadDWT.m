function [dwells,sampling,offsets,model] = loadDWT(dwtfname)
% LOADDWT  Load dwell-times from a QuB-format .dwt file
%     
%   [DWELLS,SAMPLING,OFFSETS,MODEL] = LoadDWT( FILENAME )
%   Extracts dwell-time information from the QuB-format file in FILENAME.
%   
%   DWELLS is a cell array with each element represents an idealized trace
%   with each row having the class number and dwell-time of a dwell.
%   Class numbers are 1-based and times are in frames.
%   
%   SAMPLING is the time resolution milliseconds.
%  
%   OFFSETS is a vector of the start of each trace in milliseconds relative
%   to the beginning of the file, as if it were one long trace.
%
%   MODEL is a cell array with each element having a matrix in which 
%   each row is the mean FRET values and noise stdev for each state class.
%
%   NOTE: empty idealizations (with no dwells) will not be loaded.
%
%   See also: saveDWT.

%   Copyright 2007-2017 Cornell University All Rights Reserved.



% Prompt the user for the input file if none given.
if nargin<1 || isempty(dwtfname),
    [f,p] = uigetfile( {'*.dwt','Dwell-Time Files (*.dwt)'; ...
                        '*.*','All Files (*.*)'}, 'Load a dwell-time file');
    if isequal(f,0), return; end  %user hit cancel
    dwtfname = fullfile(p,f);
end


% Find segement header lines
text = fileread(dwtfname);
[first,last,match] = regexp(text,'Segment[^\n]*','start','end','match');
first(end+1) = numel(text)+1;  %special case last dwell
nTraces = numel(match);


% Parse segment header lines
data = textscan( strjoin(match,'|'), ...
    'Segment: %f Dwells: %f Sampling(ms): %f Start(ms): %f ClassCount: %f %[^|]', ...
    'delimiter','|' );
[segid,~,sampling,offsets,~,modelText] = data{:};

assert( all(sampling==sampling(1)), 'Imaging time resolution not consistent' );
sampling = sampling(1);
offsets = offsets'/sampling;


% Parse model parameters in header lines (mean FRET and noise stdev).
nModels = numel(unique(modelText));
if nModels==1,
    % If all are the same, return just one. Many functions assume this.
    model = parseModel(modelText{1});
else
    model = cell(nTraces,1);
    for s=1:nModels,
        model{segid(s)} = parseModel(modelText{s});
    end
end


% Parse dwell-times, using 1-based class numbering and frames for time.
dwells = cell(1,nTraces);

for s=1:nTraces,
    block = text(last(s)+1:first(s+1)-1);
    data = sscanf( block, '%f' );
    dwells{segid(s)} = [data(1:2:end)+1 data(2:2:end)/sampling];
end


% Remove empty dwells.
sel = ~cellfun(@isempty,dwells);
dwells = dwells(sel);
offsets = offsets(sel);
if nModels>1, model=model(sel); end


end



function m = parseModel(text)
% Parse and shape the list of FRET parameters.
m = sscanf(text,'%f');
m = [m(1:2:end) m(2:2:end)];
end




