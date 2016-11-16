function resultTree = qub_milOptimize( dwtFilenames, modelFilename, options )
% qub_milOptimize   Optimize kinetic model using QuB's MIL algorithm
%
%   RESULT = qub_milOptimize( FILE, MODEL, options )
%   returns a QuB_Tree structure (RESULT) containing the result of using
%   the maximum likelihood model optimization algorithm (MIL) implemented
%   in QuB (qub.buffalo.edu). Dwell-time data in FILE is used by MIL to
%   optimize the starting MODEL until convergence.
%
%   OPTIONS is a structure that can include the following fields:
%    * maxIter: Maximum number of interations (int)

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% TODO:
%  - allow user to specify MIL parameter values
% Changes needed to miltreeiface:
%  - standardize and document
%  - fflush so that intermediate results are output
%  - figure out why zeros are being returned...


% Process input arguments
if nargin<2,
    error('Not enough input arguments');
end
if nargin<3,
    options = struct([]);
end

% Verify input settings
if isfield(options,'use_segments') && strcmp(options.use_segments,'together'),
    error('MIL Together is only supported method.');
end



% Load initial model
model = qub_loadTree( modelFilename );
model.data = 'Initial';


%----- Load configuration settings
if isfield(options,'maxIter')
    maxIter = options.maxIter;
else
    maxIter = 100;
end

config.ModelFile = model;
config.MaxIterations.data = int32(maxIter);
config.GroupViterbi.data  = int32(0);
config.join_segments.data = int32(0);
config.ChannelIndex.data  = int32(0);
config.ThreadCount.data   = int32(1);
config.use_segments.data  = 'together';
config.SearchLimit.data   = 10.0;
config.DFP.MaxStep.data        = 0.05;
config.DFP.ConvLL.data         = 0.0001;
config.DFP.ConvGrad.data       = 0.01;
config.DFP.MaxIterations.data  = int32(maxIter);
config.DFP.MaxRestarts.data    = 0;

if maxIter>0,
    config.Mode.data = 'optimize';
else
    config.Mode.data = 'evaluate';
end

qub_saveTree( config,'.milconfig.qtr','Properties' );



%----- Run the external MIL interface program

% Find the MIL interface program, verify there's only one
if isunix,
    milFilename = locate('miltreeiface');
elseif ispc
    milFilename = locate('miltreeiface.exe');
elseif ismac
    error('Macs are not yet supported.');
end

if numel(milFilename)<1,
    error('miltreeiface program not found. Make sure it is in your path.');
elseif numel(milFilename)>1,
    warning('More than one miltreeiface program on path.')
    milFilename = milFilename{1};
end

% Compile a job queue of MIL commands to run
if ~iscell(dwtFilenames), dwtFilenames={dwtFilenames}; end
nDwtFiles = numel(dwtFilenames);

commands    = cell(nDwtFiles,1);
outputFiles = cell(nDwtFiles,1);

for i=1:nDwtFiles,
    % Save idealization data for processing by this instance.
    outputFiles{i} = sprintf('.milresult%d.qtr',i);

    % Setup and run the MIL interface shell command
    commands{i} = sprintf('"%s" .milconfig.qtr "%s" "%s"', ...
                          milFilename, dwtFilenames{i}, outputFiles{i} );

    % For UNIX/Linux systems, the locations of the supporting QuB
    % shared libraries must be explicitly specified, even though
    % they are in the same location as the executable. Windows
    % does not have this problem.
    if isunix
        milPath = fileparts(milFilename);
        commands{i} = ['LD_LIBRARY_PATH="' milPath '" ' commands{i}];
    end
end

% Run the job queue
result = jobQueue( commands, outputFiles );

if ~all(result),
    % At least one result failed. Remember that with the way things are
    % written now, a result that fails without saving a file will not be
    % correctly recognized - the command above will never complete!!
    error('MIL failed');
end

% Process the results
for i=1:nDwtFiles,
    resultTree(i) = qub_loadTree( outputFiles{i} );

    errorCode = resultTree(i).ErrorCode.data;
    LL = resultTree(i).LL.data;
    iter = resultTree(i).Iterations.data;
    fprintf(' * E=%d LL=%f I=%d\n',errorCode,LL,iter );

    switch errorCode
        case 0
            %success
        case -3
            warning('MIL: Exceeded maximum number of iterations.');
        otherwise
            error('MIL failed: %d',errorCode);
    end
    
end % for each DWT file.

