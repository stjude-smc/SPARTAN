function [optModel,LL] = milOptimize(dwt, model)
% Optimize rate constants in current model using QuB's maximum interval
% likelihood (MIL) algorithm.
%

optModel = copy(model);

% Clear old intermediate files
delete('mil_result.qtr','result.qmf','result.qrf');


% Remove zero-state dwells that confuse MIL.
keep = ones(numel(dwt),1);
for j=1:numel(dwt),
    states = dwt{j}(:,1);
    times  = dwt{j}(:,2);

    % Remove last dwell if in dark state. These dwells result from the
    % photobleached state, which is not considered in kinetic analysis.
    if numel(states)>0 && states(end)==1,
        states = states(1:end-1);
        times  = times(1:end-1);
    end

%     % Remove last dwell, which is cut short due to photobleaching.
%     % This prevents bias in kinetic parameter estimation because
%     % the last dwell is frequently an artificial "step" from
%     % high FRET to the dark state.
%     if numel(states)<1,
%         keep(j) = 0;
%     else
%         if times(end)<=1
%             states = states(1:end-1);
%             times  = times(1:end-1);
%         end
%     end

    % Save changes, marking empty traces for removal (keep=0)
    if numel(states)==0,
        keep(j) = 0;
    else
        dwt{j} = [states times];
    end
end
dwt = dwt( logical(keep) );


% Eventually, we want to use the Baum-Welch optimized model
% as a starting point for QuB... FIXME...
% Save starting model to temporary location...
% mfname = [tempname '.qmf'];

avgRates = cell(nFiles,1);
avgStdRates = cell(nFiles,1);

for i=1:nFiles,
    [avgRates{i},avgStdRates{i},milResults(i)] = bootstrapMIL( ...
                dwt, mfname, options.bootstrapN );
end

% Save results
resultTree(1).milResults = milResults;


%----- STEP 3: Compile results into a qubtree
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

optModel.rates = rates;

close(h);


end




%%
function [avgRates,stdRates,firstResult] = bootstrapMIL(dwt, model, sampling, nBootstrap)

if nargin<3, nBootstrap=1; end


% Load the dwell-time data and the model.
nTraces = numel(dwt);

%??
nStates = model.nStates;
X = eye(nStates);
X(1,:) = 1; X(:,1) = 1;
idx_nonzero = ~logical(X);


% Construct a set of bootstrap samples of the dwell-time data.
% These samples are then processed by MIL in parallel using the
% jobQueue function.
tempDwtNames = cell(nBootstrap,1);
fret_model = [to_col(model.mu) to_col(model.sigma)];
offsets = (0:nTraces)*10000;  %doesn't matter

for i=1:nBootstrap,
    % The first bootstrap sample is all data (no sampling) - this insures
    % that the "mean" result is always the same, as expected.
    if i == 1,
        idxBootstrap = 1:nTraces;
    else
        idxBootstrap = floor(rand(nTraces,1)*nTraces)+1;
    end
    
    tempDwtNames{i} = [tempname '.dwt'];
    saveDWT( tempDwtNames{i}, dwt(idxBootstrap), offsets, fret_model, sampling );
end

% Save model to file
mfname = 'bwmodel.qmf';
delete(mfname);
qub_saveTree( model.qubTree, mfname, 'ModelFile' );

% Run MIL
result = qub_milOptimize(tempDwtNames, mfname);

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

