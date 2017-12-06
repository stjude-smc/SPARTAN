function output = updateSpartan()
% updateSpartan  Check for new versions of SPARTAN online.
%
%   [STATUS] = updateSpartan() checks for new versions of SPARTAN from the web.
%   If there are, the user will be prompted to download the new version.
%   STATUS is true if an update is available and false otherwise.
%   Checks are made at most once per day.

%   Copyright 2007-2017 Cornell University All Rights Reserved. 

DELAY_SHORT = 1;  %days between normal checks
DELAY_LONG  = 7;  %days before reminding about known updates

COMPILED_ADDR = 'https://www.dropbox.com/sh/lkmst6vrf0nubn8/AAD9LcPqzE5cjr3se4i_x0i5a?dl=0';
SOURCE_ADDR   = 'https://www.dropbox.com/sh/wlcj8maxvp5pr8f/AADtLYcojAuyE7MJVmKT_q2la?dl=0';


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
verstr = cascadeConstants('version');


%% Check online for the latest version.
fprintf('Checking for updates to SPARTAN... ');
try
    latestVerString = strtrim( urlread('https://www.dropbox.com/s/bculsb8z6j130kg/SPARTAN_version.txt?dl=1','Timeout',10) );
    latest = versionEncode(latestVerString);
catch
    fprintf('Failed. Check the address below instead:\n');
    fprintf('http://www.scottcblanchardlab.com/software\n\n');
    return;
end

needsupdate = latest > versionEncode(verstr);
if nargout>0, output = needsupdate; end
if ~needsupdate,
    fprintf('Up to date. %s >= %s\n\n',verstr,latestVerString);
    return;
else
    fprintf('New version available: %s\n\n',latestVerString);
end


%% Ask the user if they want to download now
a = questdlg( sprintf('New version available (%s). Update now?',latestVerString), ...
              'SPARTAN update', 'Yes','Not now','Stop asking', 'Yes');
switch a
    case 'Yes'
        disp('A browser window will open with new SPARTAN version.');
        disp('If this does not happen, use the address below instead:');
        if isdeployed
            % compiled versions
            web(COMPILED_ADDR,'-browser');
            disp(COMPILED_ADDR);
        else
            % source code versions
            web(SOURCE_ADDR,'-browser');
            disp(SOURCE_ADDR);
        end
        fprintf('\nFor more information, go to the following address:\n');
        fprintf('http://www.scottcblanchardlab.com/software\n\n');
        checktime = now + DELAY_LONG;

    case 'Stop asking'
        checktime = now + Inf;
    otherwise
        % User closed the window or hit "not now"
        checktime = now + DELAY_LONG;
end


end




