function keep = fret_model_filter(data, idl, minSepIn)
% Select traces with state mean FRET close to the ensemble average
% For now, assume input is from batchKinetics...

%   Copyright 2022 Cornell University All Rights Reserved.

persistent minSep;
if isempty(minSep), minSep=0.12; end

narginchk(2,3);
assert( isa(data,'Traces'), 'First argument must be Traces object' );
assert( isnumeric(idl), 'Second argument must be idealization matrix' );
% assert( isa(data,'QubModel'), 'Third argument must be QubModel object' );

% Parameters
% FIXME: these should not be hard-coded
skipFrames = 0;     %number of frames in the beginning to ignore
useFrames  = 1000;  %number of frames to use when making histograms

keep = true( data.nTraces, 1 );

if nargin<3
    result = inputdlg( 'Minimum state separation:', 'SPARTAN', 1, {num2str(minSep)} );
    if isempty(result), return; end  %user hit cancel
    result = str2double(result{1});
    if isnan(result)
        errordlg('Invalid input value');
        return;
    else
        minSep = result;
    end
else
    minSep = minSepIn;
end

nStates = max(idl(:));
trace_means = nan(data.nTraces, nStates);

% Truncate data to select region without bleaching to avoid bias
fret = data.fret( :, skipFrames+(1:useFrames) );
idl  = idl( :, skipFrames+(1:useFrames) );

% Get state mean FRET values for every trace, ignoring zero state.
for s=1:nStates
    for i=1:data.nTraces
        segment = fret(i, idl(i,:)==s );
        if ~isempty(segment)
            trace_means(i,s) = mean( segment );
        end
    end
end

% Keep traces in which mean FRET value of a state is within distance
% parameter from mean FRET value of lower state.
% If state is not occupied, it is ignored.
deviations = abs( trace_means-nanmedian(trace_means) );  %deviation from ensemble average
% deviations = abs( trace_means - to_col(model.mu) );  %deviation from model
keep = all( deviations<=minSep|isnan(deviations), 2 );

fprintf('Selected %d of %d traces (%.0f%%) \nState mean FRET:', ...
        sum(keep), numel(keep), 100*sum(keep)/numel(keep) );
disp( nanmedian(trace_means) )

% TESTING
% figure;
% for i=1:nStates
%     subplot(1,nStates+1,i);
%     hist( trace_means(:,i), 30 );
%     title( sprintf('State %d',i) );
%     xlim([0 1]);
% end
% 
% deviations = abs( trace_means-nanmedian(trace_means) );  %deviation from ensemble average
% subplot(1,nStates+1,nStates+1);
% hist( deviations(:), 30 );
% title('Deviation');
% xlim([0 0.5]);
    
end




