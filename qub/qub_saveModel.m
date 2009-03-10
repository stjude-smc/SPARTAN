function model = qub_saveModel(modelFilename)
% qub_loadModel  Loads a model file created by QuB
%     
%   [DMODEL] = qub_saveModel( filename )
%   Saves a qub model file (.qmf)
%   See qub_createModel, qub_loadModel, qub_milOptimize, qub_skmIdealize
%
%  http://www.qub.buffalo.edu


% MEX CODE
if ~exist('qub_saveModel_MEX')
    error('Can''t find QuB compatibility layer files');
end
