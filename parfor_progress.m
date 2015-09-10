function parfor_progress(N_init,title)
% This function is designed to monitor the progress of parfor exactly like
% having a normal waitbar. But because parfor workers can't use the display
% or directly communicate, a timer is launched to monitor the progress and
% update a waitbar within the GUI process. Each time the timer executes, it
% looks for a 'parfor_progress.txt' file and updates the waitbar according
% to the contents. Each time an iteration of parfor is completed by any
% worker, it will updated the 'parfor_progress.txt' file to show that one
% more iteration has been completed.
%
% Example:
%
%    N=100;    %total number of parfor iterations
%    parfor_progress(100,'Please wait...'); %create the progress bar
%    parfor i=1:N
%        sleep(rand/2);    %computation
%        parfor_progress;  %update progress (one iteration)
%    end
%    parfor_progress(0);   %close waitbar and clean up
%
% The title can be updated during the loop by passing a string. Note that this
% will not increase the progress level!
%
%    parfor_progress( 'Starting the next step...' );
%
% NOTE: If there are many iterations that are fairly short, calling this
% function for every iteration may add significantly to execution time. To
% avoid this, consider calling on every 10th iteration, etc.
% Just be sure the number of steps in initialization is also reduced!
%
% Inspired by the parfor_progress script made by Jeremy Scheff:
% http://www.mathworks.com/matlabcentral/fileexchange/32101

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% This is the file used for interprocess communication (IPC).
% WARNING: WORKS FOR LOCAL EXECUTION ONLY! (workers must agree on tempdir)
% If using a cluster use the current directory or a place all workers and the
% GUI thread can agree on.
filename = [tempdir 'parfor_progress.txt'];


% ============  Executed by worker threads only  ============ %
% If the function is called with no arguments, this is a worker telling
% us one parfor iteration has completed.
if nargin==0,
%     if ~exist(filename,'file'),
%         error('parfor_progress was not initialized properly or threads do not agree on IPC file location');
%     end
    
    % Append a 1 to the file, indicating one parfor iteration has completed.
    f = fopen(filename, 'a');
    fprintf(f, '1\n');
    fclose(f);
    return;
end


% ============  Executed by GUI thread only  ============ %
% These are only persistent within the GUI thread.
% Valid only because there can only be one parfor loop running at a time.
persistent timer_handle;  
persistent wbh;

% If just text is given, update the waitbar title.
if nargin==1 && ischar(N_init) && ~isempty(timer_handle) ...
             && isvalid(timer_handle) && ishandle(wbh),
    timer_update(timer_handle,N_init);
    return;
end
    

% Whether we are closing the progress bar or opening a new one,
% need to close and clean up the previous progress bar.
if ~isempty(timer_handle) && isvalid(timer_handle),
    stop(timer_handle);
    pause(0.01); %wait for callback to finish
    delete(timer_handle);
end
if ~isempty(wbh) && ishandle(wbh),
    close(wbh);
end
if exist(filename,'file'),
    delete(filename);
end
    
% User asked to close the progress bar. Since it was already closed above, exit.
if N_init==0,
    return;
end


% User asked to create a new progress bar.
assert( isnumeric(N_init) & N_init>0, 'Invalid initialization' );
if nargin<2 || isempty(title),
    title='Please wait...';
end
wbh = waitbar( 0, title );

% Initialize the thread communication file.
f = fopen(filename, 'w');
fprintf(f, '%d\n', ceil(N_init)); % Save N at the top of file
fclose(f);

% Start a timer that will update the progress bar every second, reading the
% IPC file to determine how many iterations have been completed.
timer_handle = timer( 'ExecutionMode','fixedSpacing', 'Period',1.0, ...
                      'BusyMode','drop', 'Name','parfor_progress', ...
                      'TimerFcn',@timer_update, 'UserData',wbh );
start(timer_handle);



end %function parfor_process



function timer_update(timer_handle,title)
% This function is called at regular intervals by a timer. Each time, it
% looks for more progress indicators written to the 'parfor_progress.txt'
% file and updates the waitbar to reflect how many iterations of parfor
% have been completed so far.
%

wbh = get(timer_handle,'UserData');
filename = fullfile(tempdir,'parfor_progress.txt');

% Verify that the IPC file exists. If not, something probably went wrong
% and the timer has become a zombie and should be killed.
if ~exist(filename,'file'),
    warning('Missing parfor_progress IPC file. Stopping timer.');
    if ~isempty(wbh) && ishandle(wbh),
        close(wbh);
    end
    stop(timer_handle);
    delete(timer_handle);
    return;
end

% Verify waitbar hasn't been closed. Usually happens when something crashed.
if ~ishandle(wbh),
    %warning('Parfor_progress waitbar closed unexpectedly. Stopping timer.');
    stop(timer_handle);
    delete(timer_handle);
    return;
end

% Load IPC file and count the fraction of iterations that have completed.
f = fopen( filename, 'r' );
progress = fscanf(f, '%d');
fclose(f);
percent = (length(progress)-1)/progress(1);
percent = max(0, min(1,percent) );

% Update the waitbar.
if nargin>1 && ischar(title),
    waitbar( percent, wbh, title );
else
    waitbar( percent, wbh );
end

end

