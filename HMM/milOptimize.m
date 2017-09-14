function [optModel,LL] = milOptimize(dwtfname, model, options)
% Optimize rate constants in current model using QuB's maximum interval
% likelihood (MIL) algorithm.

% FIXME: can probably be merged with qub_milOptimize.

narginchk(2,3);

% Remove photobleached (last zero-state) dwell, which confuses MIL.
% FIXME: this should be an option or somehow be more user controllable.
[dwt,sampling,offsets,fret_model] = loadDWT(dwtfname);
for j=1:numel(dwt),
    if size(dwt{j},1)>1 && dwt{j}(end,1)==1,
        dwt{j} = dwt{j}(1:end-1,:);
    end
end

% Save all inputs for MIL to file, excluding empty traces.
tempDwtName = [tempname '.dwt'];
sel = ~cellfun(@isempty,dwt);
saveDWT( tempDwtName, dwt(sel), offsets(sel), fret_model, sampling );

mfname = 'bwmodel.qmf';
delete(mfname);
qub_saveTree( model.qubTree, mfname, 'ModelFile' );

% Run MIL
result = qub_milOptimize(tempDwtName, mfname, options);
optModel = qub_loadModel(result.ModelFile);
LL = result.LL.data;

% % Check validity of output (non-zero FRET state rates being zero
% X = eye(model.nStates);
% X(1,:) = 1; X(:,1) = 1;
% idx_nonzero = ~logical(X);
% 
% if any( model.rates(idx_nonzero) < 10e-9 )
%     warning('Key rate estimated as zero.');
% end
% if any( model.rates(~logical(eye(model.nStates))) > 2*1000/data.sampling )
%     warning('Rate estimate way out of range.');
% end

% Delete temporary files
delete(tempDwtName);
delete('.milresult*');

% Round all rates to four significant digits for display
optModel.rates = round(optModel.rates,4,'significant');

end

