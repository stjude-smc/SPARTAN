function [dwtFilenames,skmModels,LLout] = runParamOptimizer(model,dataFilenames,options)
% batchKinetics: run parameter optimization

nFiles = numel(dataFilenames);
dwtFilenames = cell( size(dataFilenames) );

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


%----- STEP 1: Optimize params using SKM and idealize
wbh = waitbar( 0, ['Running ' options.idealizeMethod '...'] );
LLout = zeros(nFiles,1);

for i=1:nFiles
    % Load data
    data = loadTraces(dataFilenames{i});
    sprintf('%d: %s', i,dataFilenames{i});

    % Idealize data using user-specified algorithm...
    switch upper(options.idealizeMethod)
    case upper('Segmental k-means'),
        [dwt,optModel,LL,offsets] = skm( data.fret, data.sampling, model, skmOptions );

    case upper('Baum-Welch'),
        skmOptions.seperately = false;  %individual fitting not supported yet.
        [optModel,LL] = BWoptimize( data.fret, data.sampling, model, skmOptions );

        % Idealize using optimized parameters
        fretModel = [to_col(result.mu) to_col(result.sigma)];
        [dwt,~,offsets] = idealize( data.fret, fretModel, result.p0, result.A );

    case upper('ebFRET'),
        [idl,optModel,LL] = runEbFret(data.fret, model);            
        [dwt,offsets] = idlToDwt(idl);

    case upper('Thresholding'),
        error('Thresholding not implemented')
        %[dwt,offsets] = tIdealize(data, model);
        %optModel = model;
        %LL=0;
        
    case upper('MIL'),
        [optModel,LL] = milOptimize(dwt, model);

    otherwise
        error('Analysis method "%s" not recognized',options.idealizeMethod);
    end

    % Average ensemble parameters to get an average model output.
    if skmOptions.seperately,
        skmModels(i).mu    = mean( [optModel.mu], 2 );
        skmModels(i).sigma = mean( [optModel.sigma], 2 );
        LLout(i) = mean(LL);
    else
        assert( numel(optModel)==1 );
        skmModels(i) = optModel;
        LLout(i) = LL;
    end

    % Save the idealization.
    if ~strcmpi(options.idealizeMethod,'MIL')
        [p,n] = fileparts(dataFilenames{i});
        dwtFilenames{i} = fullfile( p, [n '.qub.dwt'] );
        fretModel = [to_col(skmModels(i).mu) to_col(skmModels(i).sigma)];
        saveDWT( dwtFilenames{i}, dwt, offsets, fretModel, data.sampling );
    end

    waitbar(i/nFiles,wbh);
    drawnow;
end

close(wbh);

end  %FUNCTION runParamOptimizer



