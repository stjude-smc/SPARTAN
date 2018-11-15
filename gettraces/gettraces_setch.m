function settings = gettraces_setch(settings, fieldID, input)
% Alter imaging profile settings
%
%   SETTINGS = gettraces_setch(SETTINGS,IDX,INPUT)
%   Update image profile struct SETTINGS with new parameter values given
%   in the struct INPUT. IDX is the linear index into SETTINGS.geometry,
%   which defines the subfield layout, of the field to change.
%
%   See also: gettraces_gui, cascadeConstants.

%   Copyright 2016-2018 Cornell University All Rights Reserved.

narginchk(3,3);
nargoutchk(1,1);
assert( issorted(settings.wavelengths) );


% Linear index into params.geometry for each channel
[val,idx] = sort( settings.geometry(:) );
idxFields = to_row( idx(val>0) );

% List of channel properties to be updated (crosstalk is a special case)
fnames = {'scaleFluor','chNames','chDesc','wavelengths'};

% Get index of channel associated with selected field, if any.
idxCh = settings.geometry(fieldID);


% Remove any channel associated with the selected field
if idxCh>0  %isempty(input) || ~isempty(idxCh)
    
    % Remove all entries in crosstalk matrix associated with selected channel.
    if size(settings.crosstalk,1)<=2,
        settings.crosstalk = [];
    else
        settings.crosstalk(idxCh,:) = [];
        settings.crosstalk(:,idxCh) = [];
    end
    
    % Remove linear properties associated with selected channel
    for i=1:numel(fnames)
        settings.(fnames{i})(idxCh) = [];
    end
    idxFields(idxCh) = [];
end


% (Re-)insert channel into settings struct in wavelength order.
if ~isempty(input)
    
    % Determine insertion position
    idxCh = find(input.wavelengths<settings.wavelengths, 1,'first');
    if isempty(idxCh), idxCh=numel(settings.wavelengths)+1; end  %append
    
    % Insert new value in each property list
    for i=1:numel(fnames)
        f = fnames{i};
        settings.(f) = [settings.(f)(1:idxCh-1) input.(f) settings.(f)(idxCh:end)];
    end
    idxFields = [idxFields(1:idxCh-1) fieldID idxFields(idxCh:end)];
    
    if isempty(settings.crosstalk),
        settings.crosstalk = zeros(2);
    else
        nc = size(settings.crosstalk);
        settings.crosstalk = [settings.crosstalk(1:idxCh-1,:); zeros(1,nc(2));   settings.crosstalk(idxCh:end,:)];
        settings.crosstalk = [settings.crosstalk(:,1:idxCh-1)  zeros(nc(1)+1,1)  settings.crosstalk(:,idxCh:end)];
    end
end

% Reconstruct field geometry channel assignment matrix
settings.geometry(:) = 0;
settings.geometry(idxFields) = 1:numel(idxFields);


end %function gettraces_setch



