function [ok,msg] = qub_verifyModel(model)
% qub_createModel  Creates a simple starting point model
%     
%   [MODEL] = qub_verifyModel( MODEL )
%   Checks correctness of MODEL structure
%   See also qub_saveModel, qub_loadModel, qub_milOptimize, qub_skmIdealize
%
%
%  http://www.qub.buffalo.edu

ok = false;

if ~all( isfield(model,{'mu','sigma','rates','p0'}) )
    msg = 'Not all model parameters provided - consider using qub_createModel()';
    return;
end

nStates = numel(model.mu);

if ~( nStates==numel(model.sigma)  && ...
      nStates==numel(model.p0)     && ...
      nStates==size(model.rates,1) && ...
      nStates==size(model.rates,1) ),
    msg = 'Size mismatch';
    return;
end

if abs(sum(model.p0)-1) > 0.01
    msg = 'p0 values not normalized';
    return;
end

if isfield(model,'pt')

    if numel(model.pt)~=nStates,
        msg = 'Size mismatch (pt)';
        return;
    elseif abs(sum(model.pt)-1) > 0.01,
        msg = 'pt values not normalized';
        return;
    end
    
end

if any( model.rates<0 ),
    msg = 'Negative rates not allowed.';
    return;
end
    
msg = 'Okay.';
ok = true;

