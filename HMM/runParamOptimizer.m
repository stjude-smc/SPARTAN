function [dwtfile,outModel,LL] = runParamOptimizer(model,trcfile,options)
% batchKinetics: run parameter optimization

% Setup algorithm settings
skmOptions.maxItr = options.maxItr;
skmOptions.convLL = 1e-4;
skmOptions.seperately = options.seperately;
if skmOptions.seperately,
    skmOptions.quiet = 1;  %skm  is faster than launching the waitbar
end

% Remove intermediate files from previous runs.
warning('off','MATLAB:DELETE:FileNotFound');
delete('resultTree.mat','bwmodel.qmf');

% Load data
data = loadTraces(trcfile);

if isfield(options,'dataField') && isfield(data, options.dataField)
    input = data.(options.dataField);
    
    % Normalize fluorescence intensities to fall in ~[0,1].
    % FIXME: will not work for transient events.
    % Additional options may be needed to control this behavior.
    if isempty(strfind(options.dataField,'fret'))
        temp = input(:,1:10);
        input = input / mean(temp(:));
    end
else
    input = data.fret;
end

% Idealize data using user-specified algorithm...
switch upper(options.idealizeMethod)
case upper('Segmental k-means'),
    [dwt,optModel,LL,offsets] = skm( input, data.sampling, model, skmOptions );

case upper('Baum-Welch'),
    skmOptions.seperately = false;  %individual fitting not supported yet.
    [optModel,LL] = BWoptimize( input, data.sampling, model, skmOptions );

    % Idealize using optimized parameters
    fretModel = [to_col(optModel.mu) to_col(optModel.sigma)];
    [dwt,~,offsets] = idealize( input, fretModel, optModel.p0, ...
                                        optModel.calcA(data.sampling) );

case upper('ebFRET'),
    [idl,optModel,LL] = runEbFret(input, model);            
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
    
    %outModel.mu = mean( [optModel.mu], 2 );   
    temp = arrayfun(@(x)to_col(x.mu), optModel, 'Uniform',false);
    outModel.mu = mean( cat(2,temp{:}), 2 );
    
    %outModel.sigma = mean( [optModel.sigma], 2 );
    temp = arrayfun(@(x)to_col(x.sigma), optModel, 'Uniform',false);
    outModel.sigma = mean( cat(2,temp{:}), 2 );
    
    outModel.rates = mean( cat(3,optModel.rates), 3 );
    LL = mean(LL);
else
    assert( numel(optModel)==1 );
    outModel = optModel;
end

% Truncate values to rates to four significant figures for display
outModel.rates = round(outModel.rates,4,'significant');

% Save the idealization, deleting any previous ones.
[p,n] = fileparts(trcfile);
dwtfile = fullfile( p, [n '.dwt'] );
if exist(dwtfile,'file'), delete(dwtfile); end

dwtfile = fullfile( p, [n '.qub.dwt'] );
fretModel = [to_col(outModel.mu) to_col(outModel.sigma)];
saveDWT( dwtfile, dwt, offsets, fretModel, data.sampling );

    
end  %FUNCTION runParamOptimizer



