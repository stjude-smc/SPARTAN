function result = jobQueue( commands, outputFiles, timeout )
% JOBQUEUE   Execute external programs in parallel.
%
%   RESULT = JOBQUEUE( COMMANDS, OUTFILES, timeout )
%   executes a set of shell COMMANDS (cell array of strings), each of which
%   is expected to produce a corrosponding file (cell array of strings),
%   the names of which are specified in OUTFILES (cell array of filenames).
%   TIMEOUT (optional) specifies the maximum execution time - if this time
%   is exceeded, the function will return a value, but the external program
%   may still be running.

% TODO:
%  * Automatically determine number of processes to use 

% USER TUNABLE PARAMETERS:

% number of commands to run simultaneously.
if exist('maxNumCompThreads','builtin'),
    numProcesses = maxNumCompThreads;
else
    constants = cascadeConstants;
    numProcesses = constants.nProcessors;
end



% Handle input files
nCommands = numel(commands);
assert( nCommands == numel(outputFiles) );
assert( iscell(commands) && iscell(outputFiles) );

% Remove all output files before starting
for i=1:nCommands,
    delete(outputFiles{i});
end

% Create a queue of commands to run
queue = struct( 'command',commands, 'outputFile',outputFiles, ...
                'isRunning',0, 'isDone',0 );

% Until all commands have completed, ...
tic; %keep track of start time.
while( ~all( [queue.isDone] ) )
    
    % Check the status of running commands to see if any have finished.
    runningUnits = find( [queue.isRunning] );
    
    for unitID=runningUnits(:)',        
        % If output file exists, we assume the command completed,
        % meaning that a new execution unit can start.
        if exist(outputFiles{unitID},'file')
            queue(unitID).isRunning = 0;
            queue(unitID).isDone = 1;
            disp(['Job ' num2str(unitID) ' is done.']);
        end
    end
    
    % If a command has finished, launch another to fill its place
    if sum([queue.isRunning]) < numProcesses,

        idleUnits = find( ~[queue.isDone] & ~[queue.isRunning] );
        
        unitsToStart = numProcesses-sum([queue.isRunning]);
        unitsToStart = min(unitsToStart,numel(idleUnits));
        idleUnits = idleUnits(1:unitsToStart);
        
        % Loop through all commands in the queue, checking if a command has
        % finished.
        for unitID=idleUnits(:)',
            % Launch the command listed in this unit. The '&' will cause
            % the command to be run in the background.
            system( [commands{unitID} ' &'] );
            queue(unitID).isRunning = 1;
            disp(['Started job #' num2str(unitID)]);
        end

    end

    % Only check for changes every 100ms to avoid taxing the CPU.
    pause(0.1);
    
    if nargin>=3 && toc>timeout,
        warning('Execution did not finish before timeout!');
        break;
    end
end

% Return an array containing the status of each run.
result = [queue.isDone];



