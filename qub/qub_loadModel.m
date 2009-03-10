function model = qub_loadModel(modelFilename)
% qub_loadModel  Loads a model file created by QuB
%     
%   [DMODEL] = qub_loadModel( filename )
%   Loads a qub model file (.qmf) -- specifically, this converts the
%   QUB_Tree object imbedded in a .qmf file into a structure, as
%   defined in qub_createModel.
%
%  http://www.qub.buffalo.edu


% MEX CODE
if ~exist('qub_loadModel_MEX')
    error('Can''t find QuB compatibility layer files');
end
