function varargout = removeBlinks( varargin )
% REMOVEBLINKS  Remove dark state dwells from dwell-time series (.dwt files)
% 
%   removeBlinks(FILENAME) loads .dwt file FILENAME, remove all dwells in the
%   dark state, presumed to be the lowest class number, and save the modified
%   dwell-times in a new .dwt file, prompting the user for the name.
%
%   removeBlinks() prompts the user for the input dwell-time filename.
%
%   [DWT,OFFSETS] = removeBlinks(DWT,OFFSETS) removes all dark state dwells 
%   from the dwell-times in DWT, updating OFFSETS if necessary. 
%
%   See also: correctedDwelltimes, loadDwt.

%   Copyright 2015 Cornell University All Rights Reserved.


%% Process input parameters.

% If no input given, ask for a .dwt file.
if nargin<1,
    dwtfname = getFile('*.dwt','Select a dwell-time file to process');
    if isempty(dwtfname), return; end  %user hit cancel
    
    [dwt,sampling,offsets,model] = loadDWT(dwtfname);

% Filename of .dwt file is given.
elseif ischar(varargin{1}),
    dwtfname = varargin{1};
    [dwt,sampling,offsets,model] = loadDWT(dwtfname);

% Dwell-time information is given in parameters.
elseif iscell(varargin{1}),
   if nargin<2,
       error('Not enough input arguments');
   end
   
   dwtfname = '';
   [dwt,offsets] = deal( varargin{:});
end


%% Remove dark state from dwell-times.
% We assume there can never be sequences of dark-state dwells (dwells in the
% same state) and that the dark state has the lowest class number (1).
% Times from loadDWT are in frames and must be discrete.

% Build list of dwell times in each class
for i=1:numel(dwt),
    
    % Input dwell-time data for the current trace.
    classes = dwt{i}(:,1);
    times   = dwt{i}(:,2);
    nDwells = numel(times);
    
    for j=1:nDwells,
       
        % Ignore non-zero FRET states.
        if classes(j)>1, continue; end
        
        % If the first dwell is a blink, increase the offset to skip it.
        if j==1,
            offsets(i) = offsets(i) + times(j);
            
        % Remove internal blinks. Add the time to the surrounding non-zero 
        % state dwells equally. Slight bias towards foward dwell is negligable.
        elseif j<nDwells,
            times(j-1) = times(j-1) + floor(times(j)/2);
            times(j+1) = times(j+1) +  ceil(times(j)/2);
        end
        times(j) = 0;
        
    end %for each dwell
    
    
    % Rebuild dwell-times to remove empty blink dwells and merge dwells that
    % were broken up by blinks.
    if sum(times)>0
        dwt(i) = idlToDwt(  dwtToIdl([classes times])  );
    else
        dwt{i} = [];
    end
                
end %for each trace

% TODO: should be possible to use dwtToIdl/idlToDwt at the end, but this
% seems to scramble the traces...

% Remove any empty traces that only had dark state dwells.
idxRemove = cellfun( @isempty, dwt );
dwt(idxRemove) = [];
offsets(idxRemove) = [];




%% Save the resulting corrected dwell-times.

% Not output parameters given, so we assume the user wants to save the result.
if nargout==0 && ~isempty(dwtfname),
    [p,outname,ext] = fileparts(dwtfname);
    outname = fullfile( p, [outname '_noblinks' ext] );
    saveDWT(outname,dwt,offsets,model,sampling);
else
    varargout = {dwt,offsets};
end




