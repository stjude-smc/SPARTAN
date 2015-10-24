function autoplots( direct )
% AUTOPLOTS  Automaticly display ensemble plots from new files.
%
%   AUTOPLOTS(DIR) runs MAKEPLOTS with all of the "auto.traces" files in
%   the directory DIR. New files are added automatically as they appear.
%   If no directory is given, the user will be prompted for one.
%
%   See also: makeplots.

%   Copyright 2015 Cornell University All Rights Reserved.


%%

if nargin<1,
    direct = uigetdir(pwd, 'Choose directory with data to display');
    if direct==0, return; end  %user hit cancel
end

% Ask the user how plots should be ordered.
choices = {'Alphabetically','By date modified','Most recent only'};
sel = questdlg( 'How should files be ordered?','Makeplots order', ...
                choices{:}, choices{1} );
mode = find( strcmp(sel,choices) );

% Create the file list the first time
files = filelist(direct,mode);
hfig = makeplots( files );

% Create a timer to look for 
timer_handle = timer( 'ExecutionMode','fixedSpacing', 'Period',2.0, ...
                      'BusyMode','drop', 'Name',mfilename, ...
                      'TimerFcn',@tupdate, ...% 'StopFcn',@tstop, ...
                      'UserData',{hfig,direct,mode,files} );
start(timer_handle);

% Add a callback to stop the timer if the window is closed?


end %FUNCTION MAKEPLOTSDIR



function filenames = filelist(direct,mode)
% Get a list of all auto.traces files in the given directory

files = regexpdir(direct,'\.traces$',false);
filenames = {files.name};

if nargin>=2 && mode>1,  
    % Sort by date modified
    dates = [files.datenum];
    [~,idx] = sort(dates);
    filenames = filenames(idx);
    
    % Choose only the most recent.
    if mode==3,
        filenames = filenames(end);
    end
end

end %FUNCTION FILELIST



%%
function tupdate(timer_handle, ~)
% Look for any new files and remake the plots if there are any.

params = get(timer_handle,'UserData');
[hfig,direct,mode,oldFiles] = params{:};

% Stop the timer if the window has been closed.
if ~ishandle(hfig),
    stop(timer_handle);
    delete(timer_handle);
    return;
end

files = filelist(direct,mode);

% If there are any new files, re-create the plot with the new files.
% The new plot is positioned in roughly the same location.
if ~isequal(files,oldFiles),
    pos = get(hfig,'Position');
    close(hfig);
    
    hfig = makeplots(files);
    set(hfig,'Position',pos);
    
    set(timer_handle, 'UserData', {hfig,direct,mode,files});
end

end %function tupdate



