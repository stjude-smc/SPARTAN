function [idlout,outModel,rates] = runParamOptimizer(data, idl, model, options)
% batchKinetics: run parameter optimization
% This function performs any necessary pre-processing before launching one
% of the supporting optimization algorithms.

narginchk(4,4);
nargoutchk(2,3);
idlout = idl;
outModel = model;
rates = [];

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

% Truncate data
for i=1:data.nTraces
    if options.truncate(i)>2
        input(i, options.truncate(i):end) = 0;
    else
        options.exclude(i) = true;
    end
end
input = input( ~options.exclude, : );


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
    [idl,optModel] = runHMJP( data.getSubset(~options.exclude), model, options );
    
case 'MIL'
    % Load dwell-time information
    assert( ~isempty(idl), 'Traces must be idealized before running MIL');
    dwt = idlToDwt( idl(~options.exclude,:), true );
    dt = data.sampling/1000;

    % Run MIL, only updating model rates.
    % NOTE: optModel will have the .qubTree model values, which only reflect 
    % the model as originally loaded from file. FIXME.
    if strcmpi( options.idealizeMethod, 'MIL (Together)' )
        options.updateModel = true;
        optModel = milOptimize(dwt, dt, model, options);
        model.rates = optModel.rates;
    else
        result = milOptimizeSeparately(dwt, dt, model, options);
        rates = nan( size(result,1), size(result,2), data.nTraces );
        rates(:,:,~options.exclude) = result;
        ratehist(rates);
        return;
    end

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

if ~strcmpi( options.idealizeMethod(1:3), 'MIL' )
    idlout = zeros( data.nTraces, data.nFrames );
    idlout( ~options.exclude, :) = idl;
end
    
end  %FUNCTION runParamOptimizer



