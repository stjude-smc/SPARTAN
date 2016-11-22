function [dwtfile,outModel,LL] = runParamOptimizer(model,trcfile,options)
% batchKinetics: run parameter optimization

% Setup algorithm settings
skmOptions.maxItr = options.maxItr;
skmOptions.convLL = 1e-4;
skmOptions.seperately = options.seperately;
if skmOptions.seperately,
    skmOptions.quiet = 1;
end

% Remove intermediate files from previous runs.
warning('off','MATLAB:DELETE:FileNotFound');
delete('resultTree.mat','bwmodel.qmf');

% Load data
data = loadTraces(trcfile);

% Idealize data using user-specified algorithm...
switch upper(options.idealizeMethod)
case upper('Segmental k-means'),
    [dwt,optModel,LL,offsets] = skm( data.fret, data.sampling, model, skmOptions );


case upper('Baum-Welch'),
    skmOptions.seperately = false;  %individual fitting not supported yet.
    [optModel,LL] = BWoptimize( data.fret, data.sampling, model, skmOptions );

    % Idealize using optimized parameters
    fretModel = [to_col(optModel.mu) to_col(optModel.sigma)];
    [dwt,~,offsets] = idealize( data.fret, fretModel, optModel.p0, ...
                                        optModel.calcA(data.sampling) );

case upper('ebFRET'),
    [idl,optModel,LL] = runEbFret(data.fret, model);            
    [dwt,offsets] = idlToDwt(idl);

case upper('Thresholding'),
    error('Thresholding not implemented')
    %[dwt,offsets] = tIdealize(data, model);
    %optModel = model;
    %LL=0;

otherwise
    error('Analysis method "%s" not recognized',options.idealizeMethod);
end

% Average ensemble parameters to get an average model output.
if numel(optModel)>1,
    outModel = copy(model);
    outModel.mu    = mean( [optModel.mu], 2 );
    outModel.sigma = mean( [optModel.sigma], 2 );
    outModel.rates = mean( cat(3,optModel.rates), 3 );
    LL = mean(LL);
else
    assert( numel(optModel)==1 );
    outModel = optModel;
end

% Save the idealization.
[p,n] = fileparts(trcfile);
dwtfile = fullfile( p, [n '.qub.dwt'] );
fretModel = [to_col(optModel.mu) to_col(optModel.sigma)];
saveDWT( dwtfile, dwt, offsets, fretModel, data.sampling );

    
end  %FUNCTION runParamOptimizer



