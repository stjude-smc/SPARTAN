function parfor_progress(N_init,title,wbh_init)
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
%    parfor_progress(100,'Please wait...'); %initialize
%    parfor i=1:N
%        sleep(rand/2);    %computation
%        parfor_progress;  %update progress
%    end
%    parfor_progress(0);   %close waitbar and clean up
%
%
% NOTE: only works if all workers are on the current host.
%
% NOTE: If there are many iterations that are fairly short, calling this
% function for every iteration may add significantly to execution time. TO
% avoid this, consider calling on every 10th or 100th (etc) iteration.
% Just be sure the number of steps in initialization is also reduced!
%
% Created by: Daniel Terry. Copyright 2015, Weill-Cornell Medical College
% Inspired by the parfor_progress script made by Jeremy Scheff:
% http://www.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
% 

% TODO: allow user to pass a starting (minimum) progres value and a waitbar
% handle so that a waitbar used as part of a larger block of code can be
% updated within the parfor loop seamlessly.
% Look for progress bars that are not closed before starting the next one
% and give a warning?
% Also want to allow the progress bar title to be updated somehow.


% This is the file used for interprocess communication (IPC).
% Workers must agree on tempdir! This is for local execution only!
% Replace tempdir with a place the GUI and workers will agree on if trying
% to execute on a cluster. Using the current directory usually works.
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
    
% User asked to close the progress bar.
if strcmp(N_init,'close') || N_init==0,
    return;  %nothing more to do
end


% User asked to create a new progress bar.
assert( isnumeric(N_init) & N_init>0, 'Invalid initialization' );
if nargin<2 || isempty(title),
    title='Please wait...';
end
wbh = waitbar( 0, title );

% Initialize the thread communication file.
f = fopen(filename, 'w');
fprintf(f, '%d\n', N_init); % Save N at the top of file
fclose(f);

% Start a timer that will update the progress bar every second, reading the
% IPC file to determine how many iterations have been completed.
timer_handle = timer( 'ExecutionMode','fixedSpacing', 'Period',1.0, ...
                      'BusyMode','drop', 'Name','parfor_progress', ...
                      'TimerFcn',@timer_update, 'UserData',wbh );
start(timer_handle);



end %function parfor_process



function timer_update(timer_handle,~)
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
end

% Load IPC file and count the fraction of iterations that have completed.
f = fopen( filename, 'r' );
progress = fscanf(f, '%d');
fclose(f);
percent = (length(progress)-1)/progress(1);
assert( percent>=0 && percent<=1, 'Invalid progress' );

% Update the waitbar.
waitbar( percent, wbh );

end


