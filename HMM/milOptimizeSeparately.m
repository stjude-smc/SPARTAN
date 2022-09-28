function [rates,optModel] = milOptimizeSeparately( dwt, dt, model, options )
% milOptimizeSeparately: run MIL model optimizer for each trace separately.
%
%   [RATES,MODELS] = milOptimizeSeparately( DWT, DT, MODEL ) uses the 
%   maximum interval likelihood (MIL) algorithm to optimize MODEL (can be
%   a path to a .model or .qmf file, or a QubModel object) using the 
%   dwell-time information in the cell array DWT. DT is the frame interval
%   in seconds RATES is a matrix with one row per trace and columns are
%   the rate constant for low to high (first column) and high to low
%   (second column). MODELS is an array of QubModel objects for each trace.
%
%   milOptimizeSeparately( DWT, DT, MODEL ) will save the rate constants in a
%   file named 'rates.mat' in the current directory.


% Check input arguments
narginchk(3,4);
nargoutchk(0,2);

assert( iscell(dwt), 'First argument should be cell array' );
if ischar(model)
    model = QubModel(model);
end
assert( isa(model,'QubModel'), 'Second argument should be a QubModel object or path to a .model file' );
    

nTraces = numel(dwt);
options.updateModel = false;
options.verbose = false;
options.UseParallel = false;
rates = nan( size(model.rates,1), size(model.rates,2), numel(dwt) );
optModel(nTraces) = QubModel;  %allocate array of model objects

wbh = parfor_progressbar( nTraces, 'Starting parallel pool...' );
gcp;
wbh.message = 'Running...';

parfor i=1:nTraces
    if isempty(dwt{i}), continue; end
    wbh.iterate(1); %#ok<PFBNS>
    
    try
        optModel(i) = milOptimize(dwt(i), dt, model, options);
        
        % Ignore traces that fail to converge
        %if exitflag~=1 && exitflag~=2, continue; end
        
        rates(:,:,i) = optModel(i).rates;
    catch
        continue;  %failed runs will have NaN values in rates matrix.
    end
end

% if nargout<1
    save('rates.mat','rates');
% end
close(wbh);

                                        
end  %function




