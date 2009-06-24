function resultTree = qub_milOptimize( dwtFilename, modelFilename )
% qub_milOptimize  DESC
%
%   DESC
%   

% Verify input arguments
if nargin<2,
    error('Both a DWT and model file must be specified');
end

if ~exist(dwtFilename,'file')
    error('DWT file does not exist');
end

if ~exist(modelFilename,'file')
    error('Model file does not exist');
end


%% Define file locations for input and output
baseName = tempname;
configFilename = [baseName '_milconfig.qtr'];
% dataFilename   = [baseName '_mildata.dwt'];
resultFilename = [baseName '_milresult.qtr'];

% copyfile(dwtFilename,dataFilename);
dataFilename = dwtFilename;


%% Construct a MIL input configuration tree

% Load model file
settings.ModelFile = qub_loadTree( modelFilename );
settings.ModelFile.data = 'Initial';

% Setup default optimization parameters...
settings.MaxIterations.data     = 50;
settings.GroupViterbi.data      = 0;
settings.('join_segments').data = 0;
settings.ChannelIndex.data      = 0;
settings.ThreadCount.data       = 0;
settings.Mode.data              = 'optimize';
settings.('use_segments').data  = 'together';
settings.SearchLimit.data       = 10.0;
settings.DFP.MaxStep.data       = 0.1;
settings.DFP.ConvLL.data        = 0.0001;
settings.DFP.ConvGrad.data      = 0.01;
settings.DFP.MaxIterations.data = 50;
settings.DFP.MaxRestarts.data   = 0;

% Save configuration file
qub_saveTree(settings,configFilename,'Properties');

configFilename = '.milconfig.qtr';

%% Run the external MIL optimizer

constants = cascadeConstants();

binLoc = constants.binaryLocation;

if strfind(computer,'PC')
    error('No Windows support yet');
elseif strfind(computer,'GLN')
    milCmd = ['LD_LIBRARY_PATH=' binLoc ' ' binLoc ...
              'miltreeiface ' configFilename ' ' dataFilename ' ' resultFilename]
end

% Run MIL
errorCode = unix( milCmd, '-echo' );

if errorCode~=0
    switch( errorCode )
    case -1
        error('miltreeiface: incorrect args');
    case -2
        error('miltreeiface: invalid or missing DWT file');
    otherwise
         error('miltreeiface: unknown error %d',errorCode);
    end
end


%% Collect the results for output...

if ~exist(resultFilename,'file'),
    error('No results recieved from MIL');
end

resultTree = qub_loadTree( resultFilename );

errorCode = resultTree.ErrorCode.data;
LL      = resultTree.LL.data;
iter    = resultTree.Iterations.data;
sprintf(' * E=%d LL=%f I=%d', errorCode,LL,iter )

if errorCode==-3
    warning('MIL: Exceeded maximum number of iterations.');
elseif errorCode~=0
    error('MIL failed with unknown error %d',errorCode);
end







