function [dwells,sampling,offsets,modelOut] = loadDWT(dwtfilename)
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
[sampling,offsets,modelOut] = deal( [] );
dwells  = {};

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
    
    ndwells   = data{2};
    sampling  = data{3};
    startTime = data{4};
    nClasses  = data{5};
    model     = str2num(data{6}{1});
    
    offsets(segid) = startTime/sampling;  %segment start time (frames)
    
    % Load data within segment
    data = textscan(fid, '%d%d');
    
    % Save dwells in segment
    %            time (frames)  State # (zero-based)
    dwells{segid} = [data{1}+1 data{2}/sampling];
    
end
fclose(fid);

% Convert model information into the expected matrix shape.
assert( mod(numel(model),2)==0, 'Bad FRET model definition' );
modelOut = zeros( numel(model)/2, 2 );
modelOut(:,1) = model(1:2:end); %FRET values
modelOut(:,2) = model(2:2:end); %standard deviations

if ~exist('sampling','var'),
    error('Malformed DWT file');
end

end  % function LoadDWT
