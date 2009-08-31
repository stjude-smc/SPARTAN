function resultTree = qub_milOptimize( dwtFilename, modelFilename )
%
%
%
%
%

% TODO:
%  - allow user to specify MIL parameter values
% Changes needed to miltreeiface:
%  - standardize and document
%  - fflush so that intermediate results are output
%  - figure out why zeros are being returned...


% Load initial model
model = qub_loadTree( modelFilename );
model.data = 'Initial';

%----- Load configuration settings
maxIter = 100;

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


%----- Save idealization data
copyfile( dwtFilename,'.mildata.dwt' );


%----- Run the external MIL interface program

% Find the MIL interface program, verify there's only one
if isunix,
    filename = locate('miltreeiface');
elseif ispc
    filename = locate('miltreeiface.exe');
elseif ismac
    error('Macs are not yet supported.');
end

if numel(filename)<1,
    error('miltreeiface program not found. Make sure it is in your path.');
elseif numel(filename)>1,
    warning('More than one miltreeiface program on path.')
end
filename = filename{1};

% Setup and run the MIL interface shell command
cmd = [filename ' .milconfig.qtr .mildata.dwt .milresult.qtr'];

% For UNIX/Linux systems, the locations of the supporting QuB
% shared libraries must be explicitly specified, even though
% they are in the same location as the executable. Windows
% does not have this problem.
if isunix
    milPath = fileparts(filename);
    cmd = ['LD_LIBRARY_PATH=' milPath ' ' cmd];
end

% Run MIL.
[retval,stdout] = system(cmd);


%----- Process the results from MIL
disp( sprintf('Executed: %s',cmd) );
disp( '------- BEGIN OUTPUT -------' );
disp( stdout );
disp( '------- END OUTPUT -------' );

switch retval
    case -1
        error('miltreeiface: incorrect args');
    case -2
        error('miltreeiface: invalid or missing DWT file');
    case 0
        %Success: load the result for later processing
        resultTree = qub_loadTree('.milresult.qtr');
    otherwise
        error('miltreeiface unknown error (%d)',retval);
end

errorCode = resultTree.ErrorCode.data;
LL = resultTree.LL.data;
iter = resultTree.Iterations.data;
disp( sprintf(' * E=%d LL=%f I=%d',errorCode,LL,iter ));

switch errorCode
    case 0
        %success
    case -3
        warning('MIL: Exceeded maximum number of iterations.');
    otherwise
        error('MIL failed: %d',errorCode);
end
            
