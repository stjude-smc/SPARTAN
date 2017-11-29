function settings = gettraces_setch(settings, idxNewField, input)
% Alter imaging profile settings
%
%   SETTINGS = gettraces_setch(SETTINGS,IDX,INPUT)
%   Sets the values of imaging profile struct SETTINGS for the field index IDX
%   to the values given in INPUT. If INPUT is empty, the field is removed.

%   Copyright 2016 Cornell University All Rights Reserved.

narginchk(3,3);
nargoutchk(1,1);

idxFields = find(settings.geometry);
idxCh = find(idxNewField==idxFields); %index into existing channel parameter list.
fnames = {'scaleFluor','chNames','chDesc','wavelengths'};
settings.geometry(idxNewField) = ~isempty(input);


% Remove channel or remove existing channel before re-inserting below.
% (this is done instead of updating to ensure the proper wavelength order)
if isempty(input) || ~isempty(idxCh)
    if size(settings.crosstalk,1)<=2,
        settings.crosstalk = [];
    else
        settings.crosstalk(idxCh,:) = [];
        settings.crosstalk(:,idxCh) = [];
    end
    
    for i=1:numel(fnames)
        settings.(fnames{i})(idxCh) = [];
    end
end

% Insert channel
if ~isempty(input)
    % Determine position in list to insert by wavelength order.
    idxCh = find(input.wavelengths<settings.wavelengths, 1,'first');
    if isempty(idxCh), idxCh=numel(settings.wavelengths)+1; end  %append
    
    % Insert new value in each record.
    input.idxFields = idxNewField;
    for i=1:numel(fnames)
        f = fnames{i};
        settings.(f) = [settings.(f)(1:idxCh-1) input.(f) settings.(f)(idxCh:end)];
    end
    
    % FIXME: what if settings.crosstalk is scalar?
    if isempty(settings.crosstalk),
        settings.crosstalk = zeros(2);
    else
        nc = size(settings.crosstalk);
        settings.crosstalk = [settings.crosstalk(1:idxCh-1,:); zeros(1,nc(2));   settings.crosstalk(idxCh:end,:)];
        settings.crosstalk = [settings.crosstalk(:,1:idxCh-1)  zeros(nc(1)+1,1)  settings.crosstalk(:,idxCh:end)];
    end
end

end %function gettraces_setch



