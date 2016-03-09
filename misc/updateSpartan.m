function output = updateSpartan()
% updateSpartan  Check for new versions of SPARTAN online.
%
%   [STATUS] = updateSpartan() checks for new versions of SPARTAN from the web.
%   If there are, the user will be prompted to download the new version.
%   STATUS is true if an update is available and false otherwise.
%   Checks are made at most once per day.

%   Copyright 2007-2015 Cornell University All Rights Reserved. 

DELAY_SHORT = 1;  %days between normal checks
DELAY_LONG  = 7;  %days before reminding about known updates


% Check the website at most once a day. Otherwise, just return the current
% version status.
persistent needsupdate;
if isempty(needsupdate), needsupdate = false; end

persistent checktime;
if ~isempty(checktime) && now<checktime,
    % Check was performed recently. Skip for now to save time.
    if nargout>0, output = needsupdate; end
    return;
end
checktime = now + DELAY_SHORT;

% Get current version number
constants = cascadeConstants;
current = cellfun( @(s)sscanf(s,'%f'), strsplit(constants.version,'.') );


%% Check online for the latest version.
fprintf('Checking for updates to SPARTAN... ');
try
    input = urlread('https://www.dropbox.com/s/bculsb8z6j130kg/SPARTAN_version.txt?dl=1');
    latest = cellfun( @(s)sscanf(s,'%d'), strsplit(input,'.') );
    assert( isnumeric(latest) || numel(latest)==3, 'Invalid version number' );
catch
    fprintf('Failed. Check the address below instead:\n');
    disp('http://www.scottcblanchardlab.com/software');
    return;
end

needsupdate = latest(1)>current(1) || latest(2)>current(2) || latest(3)>current(3);
if nargout>0, output = needsupdate; end
if ~needsupdate,
    fprintf('Up to date.\n\n');
    return;
else
    fprintf('New version available: %d.%d.%d\n\n',latest);
end


%% Ask the user if they want to download now
a = questdlg( sprintf('New version available (%s). Update now?',strtrim(input)), ...
              'SPARTAN update', 'Yes','Not now','Stop asking', 'Yes');
switch a
    case 'Yes'
        if isdeployed
            % compiled versions
            status = web('https://www.dropbox.com/sh/lkmst6vrf0nubn8/AAD9LcPqzE5cjr3se4i_x0i5a?dl=0','-browser');
        else
            % source code versions
            status = web('https://www.dropbox.com/sh/wlcj8maxvp5pr8f/AADtLYcojAuyE7MJVmKT_q2la?dl=0','-browser');
        end
        if status~=0,
            % Dropbox failed. Try the website instead
            status = web('http://www.scottcblanchardlab.com/software','-browser');
            if status~=0,
                disp('A new version of SPARTAN is available. Check the address below:');
                disp('http://www.scottcblanchardlab.com/software');
                disp('');
            end
        end

    case 'Stop asking'
        checktime = now + Inf;
    otherwise
        % User closed the window or hit "not now"
        checktime = now + DELAY_LONG;
end


end





