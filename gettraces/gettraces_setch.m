function settings = gettraces_setch(settings,idxNewField,input)
% Alter imaging profile settings
%
%   SETTINGS = gettraces_setch(SETTINGS,IDX,INPUT)
%   Sets the values of imaging profile struct SETTINGS for the field index IDX
%   to the values given in INPUT. If INPUT is empty, the field is removed.

%   Copyright 2016 Cornell University All Rights Reserved.

narginchk(3,3);
nargoutchk(1,1);

idxCh = find(idxNewField==settings.idxFields); %index into existing channel parameter list.



%% Remove a channel
if isempty(input)
    [~,idxAcceptor] = ismember(settings.chNames, {'acceptor','acceptor2'});
    if idxAcceptor>0, settings.scaleAcceptor(idxAcceptor)=[]; end
    % FIXME: acceptor2 may need to be renamed if acceptor1 is removed...

    settings.crosstalk(idxCh,:) = [];
    settings.crosstalk(:,idxCh) = [];

    fnames = {'idxFields','chNames','chDesc','wavelengths'};
    for i=1:numel(fnames)
        settings.(fnames{i})(idxCh) = [];
    end
    
    return;
end


%% Alter a channel
input.idxFields = idxNewField;
fnames = {'idxFields','chNames','chDesc','wavelengths'};

% Alter a channel in place
if ~isempty(idxCh),
    for i=1:numel(fnames)
        f = fnames{i};
        if iscell(settings.(f))
            settings.(f){idxCh} = input.(f);
        else
            settings.(f)(idxCh) = input.(f);
        end
    end
    % FIXME: alterations of wavelength should change the order?!
    
% Insert a channel
else
    % Determine position in list to insert new field
    idxCh = find(idxNewField<settings.idxFields, 1,'last');
    if isempty(idxCh), idxCh=numel(settings.idxFields)+1; end
    
    for i=1:numel(fnames)
        f = fnames{i};
        settings.(f) = [settings.(f)(1:idxCh-1) input.(f) settings.(f)(idxCh:end)];
    end
    
    nc = size(settings.crosstalk);
    settings.crosstalk = [settings.crosstalk(1:idxCh-1,:); zeros(1,nc(2));   settings.crosstalk(idxCh:end,:)];
    settings.crosstalk = [settings.crosstalk(:,1:idxCh-1)  zeros(nc(1)+1,1)  settings.crosstalk(:,idxCh:end)];
    
    % FIXME: figure out how to add acceptor scaling here.
%     [~,i] = ismember(settings.chNames, {'acceptor','acceptor2'});
%     if i>0,  settings.scaleAcceptor(i) = [];  end
end


end %function gettraces_setch


