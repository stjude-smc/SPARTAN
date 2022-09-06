function [idl,outModel] = runParamOptimizer(data, dwtfile, model, options)
% batchKinetics: run parameter optimization

narginchk(4,4);
nargoutchk(2,2);

% Remove intermediate files from previous runs.
warning('off','MATLAB:DELETE:FileNotFound');
delete('resultTree.mat','bwmodel.qmf');

% Load data
assert( isa(data,'Traces'), 'Input must be Traces object or path to .traces file' );
if ischar(data),  data=loadTraces(trcfile);  end

if isfield(options,'dataField') && ~isempty(options.dataField)
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
switch upper(options.idealizeMethod(1:3))
case 'SEG'
    [idl,optModel] = skm( input, data.sampling, model, options );

case 'BAU'
    options.seperately = false;  %individual fitting not supported yet.
    [idl,optModel] = BWoptimize( input, data.sampling, model, options );

case 'EBF'
    [idl,optModel] = runEbFret(input, data.sampling, model, options);
    
case 'MPL'
    [idl,optModel] = mplOptimize( input, data.sampling, model, options );

case 'HMJ'
    [idl,optModel] = runHMJP( data, model, options );
    
case upper('THR'),
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
    outModel.mu    = mean( cat(3,optModel.mu),    3 );
    outModel.sigma = mean( cat(3,optModel.sigma), 3 );
    outModel.p0    = mean( cat(3,optModel.p0),    3 );
    
    rates = cat(3,optModel.rates);
    outModel.rates = mean(rates,3);
    save('rates.mat','rates');
    %ratehist(rates);
else
    assert( numel(optModel)==1 );
    outModel = optModel;
end

% Truncate values to rates to four significant figures for display
outModel.rates = round(outModel.rates,4,'significant');

% Save the idealization, deleting any previous ones.
if ~isempty(dwtfile)
    saveDWT( dwtfile, idl, outModel, data.sampling );
end

    
end  %FUNCTION runParamOptimizer



