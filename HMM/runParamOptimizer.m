function [resultTree,dwtFilenames] = runParamOptimizer( model,dataFilenames,options)
% batchKinetics: run parameter optimization

resultTree = struct([]);

% disp(options);
% disp(model);
% 
% return;


h = waitbar(0,'Initializing...');

nFiles = numel(dataFilenames);
dwtFilenames = cell( size(dataFilenames) );

% Setup algorithm settings
skmOptions.maxItr = options.maxItr;
skmOptions.convLL = 1e-4;
skmOptions.seperately = options.seperately;
if skmOptions.seperately,
    skmOptions.quiet = 1;
end

bwOptions = skmOptions;
thresholdOptions = struct([]);


% Remove intermediate files from previous runs.
warning('off','MATLAB:DELETE:FileNotFound');
delete('resultTree.mat','mil_result.qtr','result.qmf','result.qrf','bwmodel.qmf');


if ~strcmp(options.idealizeMethod,'Do Nothing'),
    %----- STEP 1: Optimize params using SKM and idealize
    waitbar( 0, h, ['Idealizing data using ' options.idealizeMethod '...'] );
    
    skmLL = zeros(nFiles,1);

    for i=1:nFiles
        filename = dataFilenames{i};
        sprintf('%d: %s', i,filename);

        % Load data
        d = loadTraces(dataFilenames{i});
        data = d.fret;
        sampling = d.sampling;

        % Idealize data using user-specified algorithm...
        if strcmpi(options.idealizeMethod,'Segmental k-means'),
            
            [dwt,optModel,LL,offsets] = skm( data, sampling, model, skmOptions );
            skmLL(i) = LL(end);
            
        elseif strcmpi(options.idealizeMethod,'Baum-Welch'),
            % FIXME: this doesn't seem to actually read the result!!!
            result = BWoptimize( data, sampling, model, bwOptions );
            fretModel = [to_col(model.mu) to_col(model.sigma)];
            optModel = model;
            % TODO: update optModel with optimized parameter values.
            
            [dwt,~,offsets,LL] = idealize( ...
                    data, fretModel, result.p0, result.A );
            skmLL(i) = mean(LL);
            
        elseif strcmpi(options.idealizeMethod,'ebFRET'),
            % Verify ebFRET installed.
            if ~exist('ebfret','file')
                errormsg('ebFRET package not found. Download and add to MATLAB path.');
                return;
            end
            
            [idl,optModel] = runEbFret(data, model);            
            [dwt,offsets] = idlToDwt(idl);
            LL = 0;
            
        elseif strcmpi(options.idealizeMethod,'Thresholding'),
            
            [dwt,offsets] = tIdealize( data, model, thresholdOptions );
            optModel = model;
            
        else
            error('Analysis method "%s" not recognized',options.idealizeMethod);
        end
        
        % Remove final zero-state dwell, if it exists.
        % These dwells confuse MIL.
        keep = ones(numel(dwt),1);
        for j=1:numel(dwt),
            states = dwt{j}(:,1);
            times  = dwt{j}(:,2);
            
%             % Remove last dwell if in dark state. These dwells result from the
%             % photobleached state, which is not considered in kinetic analysis.
%             if numel(states)>0 && states(end)==1,
%                 states = states(1:end-1);
%                 times  = times(1:end-1);
%             end
%             
%             % Remove last dwell, which is cut short due to photobleaching.
%             % This prevents bias in kinetic parameter estimation because
%             % the last dwell is frequently an artificial "step" from
%             % high FRET to the dark state.
%             if numel(states)<1,
%                 keep(j) = 0;
%             else
%                 if times(end)<=1
%                     states = states(1:end-1);
%                     times  = times(1:end-1);
%                 end
%             end
            
            % Save changes, marking empty traces for removal (keep=0)
            if numel(states)==0,
                keep(j) = 0;
            else
                dwt{j} = [states times];
            end
        end
        dwt = dwt( logical(keep) );
        offsets = offsets( logical(keep) );
        
        % Save results.
        % NOTE: there is no one result when each trace is optimized seperately.
        % An "average" model is extracted instead.
        % FIXME: always average as "seperately" setting not always valid.
        if skmOptions.seperately,
            skmModels(i).mu    = mean( [optModel.mu], 2 );
            skmModels(i).sigma = mean( [optModel.sigma], 2 );
            skmModels(i).LL    = mean( LL, 2 );
        else
            skmModels(i) = optModel(1);
        end

        % Save the idealization
        [p,n] = fileparts(filename);
        dwtFilenames{i} = fullfile( p, [n '.qub.dwt'] );
        fretModel = [to_col(skmModels(i).mu) to_col(skmModels(i).sigma)];
        saveDWT( dwtFilenames{i}, dwt, offsets, fretModel, sampling );

        waitbar(0.33*i/nFiles,h);
        drawnow;
    end

    % Save results
    resultTree(1).skmModels = skmModels;
    resultTree.skmLL = skmLL;
end

% If no further action is neccessary, exit.
if strcmp(options.kineticsMethod,'Do Nothing'),
    close(h);
    return;
end




%----- STEP 2: Refine kinetic parameter estimates using MIL
waitbar(0.33,h,'Refining kinetic model using QuB...');

% Eventually, we want to use the Baum-Welch optimized model
% as a starting point for QuB... FIXME...
% Save starting model to temporary location...
% mfname = [tempname '.qmf'];
mfname = 'bwmodel.qmf';
delete(mfname);
qub_saveTree( model.qubTree, mfname, 'ModelFile' );

avgRates = cell(nFiles,1);
avgStdRates = cell(nFiles,1);

for i=1:nFiles
    [p,n] = fileparts(dataFilenames{i});
    dwtFilename = fullfile( p, [n '.qub.dwt'] );
    disp( sprintf('%d: %s', i,dwtFilename) );
    
    % Verify the idealization has been performed
    if ~exist(dwtFilename,'file'),
        error('File not idealized. Select an idealization method.');
    end

    % Run MIL
    waitbarBounds = 0.33+0.66*[(i-1) i]/nFiles;
    [avgRates{i},avgStdRates{i},milResults(i)] = bootstrapMIL( ...
                dwtFilename, mfname, options.bootstrapN );

    waitbar(waitbarBounds(end),h);
    drawnow;
end

% Save results
resultTree(1).milResults = milResults;


%----- STEP 3: Compile results into a qubtree
waitbar(1,h,'Saving results...');


rates = [];
stdRates = [];
nStates = model.nStates;

for i=1:nFiles
%     modelTree = milResults(i).ModelFile;
%     model = qub_loadModel( modelTree );
%     Q = model.rates';
    Q    = avgRates{i}';
    Qstd = avgStdRates{i}';

    idx         = find( ~logical(eye(nStates)) );
    idx_nonzero = find( ~logical(eye(nStates-1)) );

    % 1->[2,3,4], 2->[1,3,4], 3->[1,2,4], 4->[1,2,3]
    Q_parts   = Q(idx)';
    Q_nonzero = Q(2:end,2:end);
    Q_nonzero = Q_nonzero(idx_nonzero)';

    Qstd_parts   = Qstd(idx)';
    Qstd_nonzero = Qstd(2:end,2:end);
    Qstd_nonzero = Qstd_nonzero(idx_nonzero)';

    rates    = [rates    ; Q_nonzero   ];
    stdRates = [stdRates ; Qstd_nonzero];
end
nRates = size(rates,2);

% Combine rates and errors into a single matrix
output = zeros(nFiles,nRates*2);
output(:,1:2:end) = rates;
output(:,2:2:end) = stdRates;

% Construct names for each of the rates
rateNames = cell(0,1);
for i=2:nStates,
    for j=2:nStates
        if i==j, continue; end
        rateNames{end+1} = sprintf('k%d->%d', i,j );
        rateNames{end+1} = ['d' rateNames{end}];
    end
end

% Save results to file with headers
fid = fopen('rates.txt','w');

header = {'Dataset',rateNames{:}};
fprintf(fid, '%s\t', header{:});

for i=1:nFiles,
    [~,fname] = fileparts( dataFilenames{i} );
    fprintf(fid, '\n%s',  fname );
    fprintf(fid, '\t%.4f', output(i,:));
end
fclose(fid);

resultTree.rates = rates;

close(h);


end




%%
function [avgRates,stdRates,firstResult] = bootstrapMIL( ...
                          dwtFilename, modelFilename, nBootstrap )

if nargin < 3,
    nBootstrap = 1;
end


% Load the dwell-time data and the model.
[dwt,data.sampling,offsets,fret_model] = loadDWT( dwtFilename );
nTraces = numel(dwt);

nStates = numel(fret_model)/2;
X = eye(nStates);
X(1,:) = 1; X(:,1) = 1;
idx_nonzero = ~logical(X);


% Construct a set of bootstrap samples of the dwell-time data.
% These samples are then processed by MIL in parallel using the
% jobQueue function.
% idxBootstrap = cell(nBootstrap,1);
tempDwtNames = cell(nBootstrap,1);

for i=1:nBootstrap,
    % The first bootstrap sample is all data (no sampling) - this insures
    % that the "mean" result is always the same, as expected.
    if i == 1,
        idxBootstrap = 1:nTraces;
    else
        idxBootstrap = floor(rand(nTraces,1)*nTraces)+1;
    end
    
    tempDwtNames{i} = [tempname '.dwt'];
    saveDWT( tempDwtNames{i}, dwt(idxBootstrap), ...
             offsets(idxBootstrap), fret_model, data.sampling );
end

% Run MIL
result = qub_milOptimize( tempDwtNames, modelFilename );

% Process the results.
bootstrapRates = [];

for i=1:nBootstrap,
    
    % Save the first result as "the" result; others are only informative of
    % the error in the analysis.
    if i == 1,
        firstResult = result(i);
    end
    
    % Process the results for averaging
    modelTree = result(i).ModelFile;
    model = qub_loadModel( modelTree );
    rates = model.rates(:);
    
    disp( rates(:)' );
    if any( rates(idx_nonzero) < 10e-9 )
        warning( 'Key rate estimated as zero. Ignoring result.' );
        continue;
    end

    if any( rates(~logical(eye(nStates))) > 2*1000/data.sampling )
        warning( 'Rate estimate way out of range. Ignoring.' );
        continue;
    end
    
    bootstrapRates(i,:) = rates(:);
end

% Calculate average rates
avgRates = bootstrapRates(1,:);
% avgRates = mean(bootstrapRates,1);
stdRates = std(bootstrapRates,0,1);

if nBootstrap==1,
    stdRates = zeros( size(avgRates) );
end

avgRates = reshape( avgRates, size(model.rates) );
stdRates = reshape( stdRates, size(model.rates) );

save( [dwtFilename '.ratedata.txt'], 'bootstrapRates', '-ASCII' );


% Delete temporary data to save disk space.
for i=1:numel(tempDwtNames),
    delete( tempDwtNames{i} );
end
delete('.milresult*');


end


