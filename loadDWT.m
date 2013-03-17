function [dwells,sampling,offsets,model] = loadDWT(dwtfilename)
% LOADDWT  Loads a QuB idealization dwell time file
%     
%   [DWELLS,SAMPLING,OFFSETS,FRET_MODEL] = LoadDWT( DWTFILENAME )
%   Loads a .dwt file (DWTFILENAME), returning the sequence of DWELLS
%   as a cell array (Nx1, where N=number of traces). Each element has
%   an Mx2 matrix of the state number (1-based) and the duration
%   (in frames), where M=number of dwells. SAMPLING is milliseconds.
%   If no filename is given, the user is prompted for one.
%

% Set default values for outputs in case user cancels.
[sampling,offsets,nStates] = deal( [] );
dwells  = {};
model = {};

% Ask user for file if none specified.
if nargin<1 || isempty(dwtfilename),
    [f,p] = uigetfile( {'*.dwt;*.traces','Dwell-Time Files (*.dwt)'; ...
                        '*.*','All Files (*.*)'}, 'Select a DWT file');
    if f==0, return; end
    dwtfilename = [p f];
end

% Load dwell time file
if ~exist(dwtfilename,'file')
    error('Can''t find DWT file for %s',dwtfilename);
end

fid = fopen(dwtfilename);

while 1,
    
    % Load next segment in file
    % Segment: 1 Dwells: 6 Sampling(ms): 10 Start(ms): 0 ClassCount: 4 0.01
    %                                 0.034 0.15 0.061 0.3 0.061 0.55 0.061
    data = textscan(fid, 'Segment: %f Dwells: %f Sampling(ms): %f Start(ms): %f ClassCount: %f %[^\n]');
    segid = data{1};
    if numel(segid) == 0, break; end
    
    %ndwells   = data{2};
    sampling  = data{3};
    startTime = data{4};
    %nClasses  = data{5};
    
    % Reshape model to expected dimensions (FRET first column, std second).
    m = sscanf(data{6}{1},'%f');
    modelOut = zeros( numel(m)/2, 2 );
    modelOut(:,1) = m(1:2:end); %FRET values
    modelOut(:,2) = m(2:2:end); %standard deviations
    model{segid} = modelOut;
    nStates(segid) = numel(m)/2;
    
    offsets(segid) = startTime/sampling;  %segment start time (frames)
    
    % Load data within segment
    data = textscan(fid, '%f%f');
    
    % Save dwells in segment
    dwells{segid} = [data{1}+1 data{2}/sampling];
    
end
fclose(fid);

% If all the models are identical, merge into a single model. Many
% QuB-based functions assume this.
if numel( unique(nStates) )==1,
    if all(model{1}(:)==model{end}(:))
        model = model{1};
    end
end

if ~exist('sampling','var'),
    error('Malformed DWT file');
end

end  % function LoadDWT
