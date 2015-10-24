function wbh = parfor_progress(varargin)
% PARFOR_PROGRESS  Progress bar suitable for multi-process computation.
%
%   H = PARFOR_PROGRESS(N,'message', 'property',value, ...) creates a
%   progress bar that can be updated after each iteration of a PARFOR loop
%   with N total iterations before completion. A waitbar H is used to
%   display the progress level. Optional property/value pairs can be used
%   to change the appearance and behavior of the waitbar.
%
%   PARFOR_PROGRESS(H) updates the progress bar by one iteration.
%
%   PARFOR_PROGRESS(H,X) updates the progress bar by X iterations.
%
%   PARFOR_PROGRESS(H,'message') updates the task description, but not
%   progress. This should only be used outside of the PARFOR loop.
%
%   Example:
%   --------
%     N=50;    %total number of parfor iterations
%     h = parfor_progress(N,'Please wait...');  %create the progress bar
%     parfor i=1:N,
%         pause(rand);         %computation
%         parfor_progress(h);  %update progress (one iteration)
%     end
%     delete(h);   %close progress bar
%
%   Notes:
%   ------
%   Each call incurs some overhead. If there are many short iterations,
%   consider calling every 10th or 100th iteration, passing the number of
%   iterations completed between calls as a parameter.
%
%   The UserData of H is used internally. Do no modify it.
%   A temporary file is used to record the completion of each iteration.
%   A timer monitors this file and updates the waitbar accordingly.
%   This solves the problem of inter-process communication.
%
%   Inspired by the parfor_progress script made by Jeremy Scheff:
%   http://www.mathworks.com/matlabcentral/fileexchange/32101
%
%   See also: WAITBAR, PARFOR.

%   Copyright 2015 Cornell University All Rights Reserved.


narginchk(1,Inf);
nargoutchk(0,1);


%% Create a new waitbar
if isnumeric(varargin{1}),
    
    N_init = ceil(varargin{1});
    
    % Create the waitbar, passing any additional parameters on.
    if nargin>=2 && ischar(varargin{2}),
        wbh = waitbar(0,varargin{2:end});
    else
        wbh = waitbar(0,'Please wait...');
    end
    
    % Create inter-process communication file, avoiding existing files.
    % The first line is the total number of iterations expected.
    while true,
        f = sprintf('parfor_progress%d.txt', round(rand*1000));
        filename = fullfile(tempdir, f);
        if ~exist(filename,'file'), break; end
    end
    
    f = fopen(filename, 'w');
    fprintf(f, '%d\n', N_init);
    fclose(f);

    % Start a timer to regularly read the file and update the waitbar.
    timer_handle = timer( 'ExecutionMode','fixedSpacing', 'Period',1.0, ...
                          'BusyMode','drop', 'Name','parfor_progress', ...
                          'TimerFcn',@tupdate, 'StopFcn',@tstop, ...
                          'UserData',{wbh,filename} );
    start(timer_handle);
                      
    % Save internal data in the waitbar's AppData for later calls.
    setappdata(wbh,'pp_timer_handle',timer_handle);
    setappdata(wbh,'pp_filename',filename);
end


%% Update an existing waitbar
if isa(varargin{1},'handle'),
    wbh = varargin{1};
    
    % If the internal state is not valid, the waitbar was closed early.
    if ~isvalid(wbh) || ...
       ~isappdata(wbh,'pp_timer_handle') || ~isappdata(wbh,'pp_filename'),
        return;
    end
    
    timer_handle = getappdata(wbh,'pp_timer_handle');
    filename     = getappdata(wbh,'pp_filename');
    
    if ~exist(filename,'file'),
        return;
    end
    
    % Update progress by one iteration.
    if nargin==1,
        fid = fopen(filename, 'a');
        fprintf(fid, '1\n');
        fclose(fid);
    
    % Update progress by multiple iterations.
    elseif isnumeric(varargin{2}) && isscalar(varargin{2})
        fid = fopen(filename, 'a');
        fprintf(fid, '%d\n',varargin{2});
        fclose(fid);
        
    % Update the waitbar message only. Should not be called by workers.
    elseif ischar(varargin{2})
        tupdate(timer_handle,varargin{2});
        
    else
        error('Invalid arguments');
    end
end


end %function parfor_process



function tupdate(timer_handle, title)
% Function called by the timer every second to update the waitbar with the
% current progress, read from the IPC file.

params = get(timer_handle,'UserData');
[wbh,filename] = params{:};

% If the waitbar has been closed, stop the timer and clean up.
if ~ishandle(wbh) || ~exist(filename,'file'),
    stop(timer_handle);
    return;
end

% Calculate the fraction of completed iterations from IPC file.
fid = fopen( filename, 'r' );
progress = fscanf(fid, '%d');
fclose(fid);

percent = sum(progress(2:end))/progress(1);
percent = max(0, min(1,percent) );

% Update the waitbar.
if nargin>1 && ischar(title),
    waitbar( percent, wbh, title );
else
    waitbar( percent, wbh );
end

end %function tupdate



function tstop(timer_handle, ~)
% Executed if the timer crashes or if the progress bar is closed.

params = get(timer_handle,'UserData');
[wbh,filename] = params{:};

if ishandle(wbh),
    close(wbh);
end

% Wait for any callback to finish before deleting or it may crash.
pause(0.01);
delete(timer_handle);

if exist(filename,'file'),
    delete(filename);
end
    
end %function tstop




