function [idl,L] = ebFretDwt(matfile,trcfile,nStates)
% ebFretDwt  Extract dwell-time information from ebFRET session files
%
%   ebFretDwt loads session files saved by ebFRET and extracts the
%   idealization and optimal model parameters. 
%
%   ebFretDwt(FILENAME) loads the specified ebFRET session .mat file.
%   If no outputs are requested, QuB-format dwell time files (.dwt) will be
%   saved with extensions of "_ebX.qub.dwt", where X is the number of states
%   fitted in ebFRET.
%
%   [idl,models,L] = ebFretDwt() returns the state assignment matrix (idl),
%   QubModel representing the optimal model, and lower-bound evidence (L).
%   idl and model are cell arrays, one per number of fitted states, some of
%   which may be empty. L is a row vector.
%   
%   See also: loadDWT, vbFRET_dwt.

%   Copyright 2017 Cornell University All Rights Reserved.

% FIXME: ebFRET is removing a frame (!!!)


%% Process input arguments
narginchk(0,3);
nargoutchk(0,3);

% Get filename from user
if nargin<1,
    [f,p] = uigetfile('*.mat','Select a ebFRET session file');
    if ~ischar(f), return; end  %user hit cancel
    matfile = fullfile(p,f);
end
if nargin<2
    [f,p] = uigetfile('*.traces','Select the corresponding .traces file');
    if ~ischar(f), return; end
    trcfile = fullfile(p,f);
end
if nargin<3,
    a = inputdlg('Model complexity to load');
    if isempty(a), return; end  %user hit cancel
    nStates = str2double(a{1});
end


%% Load ebFRET session including analysis results
result = load(matfile);
data = loadTraces(trcfile);

assert( numel(result.series)==data.nTraces );
assert( numel(result.series(1).time)+1==data.nFrames ); %ebFRET appears to remove a frame (??)

assert( nStates<=numel(result.analysis), 'Invalid input model complexity' );
L = result.analysis(nStates).lowerbound;  %evidence per trace

% Load idealization
idl = zeros(data.nTraces,data.nFrames);  %zero means unassigned
for t=1:data.nTraces
    % Get idealization from ebFRET results
    trace = [0; result.analysis(nStates).viterbi(t).state];
    
%     % Sort idealized states so they are in order of increasing mean FRET.
%     % Leave empty states (NaN) in the current location.
%     tracemu = idlparam(data.fret(t,:), trace, nStates);
%     [~,idx] = sort( tracemu(~isnan(tracemu)) );
%     
%     idx2 = nan(size(tracemu));
%     idx2( ~isnan(tracemu) ) = idx;
%     
%     for i=1:nStates
%         trace( trace==i ) = idx2(i);
%     end    
    
    idl(t,1:numel(trace)) = trace;
end

% % Initial estimate of ensemble-average state FRET values.
% ensemblemu = idlparam(data.fret, idl, nStates);
% 
% % Reassign idealization to nearest ensemble-average states
% % This actually avoids the need to sort above in principle?
% for t=1:data.nTraces
%     trace = idl(t,:);
%     tracemu = idlparam(data.fret(t,:), trace, nStates);
%     
%     % Merging similar states
%     for i=1:nStates
%         [~,idx] = min( abs(tracemu(i)-ensemblemu) );
%         idl(t, trace==i ) = idx;
%     end
%     
% % %     Without merging similar states:
% %     used = false(size(tracemu));
% %     states = 1:nStates;
% %     for i=1:nStates
% %         % For each observed state, assign to the nearest ensemble state
% %         % FIXME?: this is not the optimal solution (depends on order)
% %         [~,idx] = min( abs(tracemu(i)-ensemblemu(~used)) );
% %         state = states(~used);
% %         idl(t,trace==i) = state(idx);
% %         used(idx) = true;
% %     end
% end

[mu,sigma] = idlparam(data.fret, idl, nStates);
disp(mu);

% Save .dwt files
if nargout==0,
    [p,f] = fileparts(trcfile);
    dwtname = fullfile(p,[f '.qub.dwt']);

    [dwt,offsets] = idlToDwt(idl);
    saveDWT(dwtname, dwt, offsets, [mu,sigma], data.sampling);
end


end



function [mu,sigma] = idlparam(fret,idl,nStates)
% Observed mean state FRET mean and noise stdev from idealization.

[mu,sigma] = deal( nan(nStates,1) );

for class=1:nStates
    edata = fret( idl==class );
    if numel(edata)>1
        mu(class) = mean(edata);
        sigma(class) = std(edata);
    end
end

end


