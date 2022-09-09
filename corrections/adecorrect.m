function data = adecorrect(input, varargin)
% Correct for acceptor direct excitation.
%
%   ADE_CORRECT( FILES, ETA_d, ETA_a )
%   For each file path in the cell array FILES, subtract acceptor direct
%   excitation from acceptor intensities, using the extinction coefficients
%   of the donor (ETA_d) and acceptor (ETA_a) to estimated the amount of
%   this fluorescence expected.
%   Creates a new output file with extension "_ade.traces" for each input.
%
%   ADE_CORRECT( FILES, FACTOR ) as above, but uses a fixed direct
%   excitation fraction (e.g., 0.05) for corrections).
%
%   NOTE: baseline subtraction in SPARTAN implicitly removes the elevated 
%   acceptor baseline due to acceptor direct excitation in cases where the 
%   donor photobleaches first, so those traces are not corrected.



%% Process input arguments
narginchk(1,3);

if nargin==3
    % See Methods in Molecular Biology, vol. 350: Protein Folding Protocols,
    %   Chapter 8 by Ben Schuler, Note #3.
    correctionFactor = 1/(1+(varargin{1}/varargin{2}));
elseif nargin==2
    correctionFactor = varargin{1};
elseif nargin==1
    correctionFactor = str2double( inputdlg('ADE correction value:',mfilename,1,{'0.05'}) );
    if isnan(correctionFactor), return; end
end

if iscell(input)
    for i=1:numel(input)
        data = adecorrect( input{i}, correctionFactor );
    end
    return;
end


%%
    if ischar(input)
        data = loadTraces(input);
    elseif isa(input,'Traces')
        data = input;
    else
        error('Input should be path to .traces file or Traces object');
    end
    endval = nan(data.nTraces,1);
    
    % Correct only traces where acceptor photobleaches first.
    % Traces in which the donor bleaches first are handled "accidentially"
    % by SPARTAN's baseline subtract procedure.
    for n=1:data.nTraces
        trace = data.fret(n, data.fret(n,:)~=0 );
        if numel(trace)>5
            endval(n) = mean(trace(end-5:end));
        end
    end
    
    idx = endval < 0.075;
    fprintf( 'Donor bleaches first in %d of %d (%0.0f%%)\n\n', ...
             nansum(idx), sum(~isnan(idx)), 100*nansum(idx)/sum(~isnan(idx)) );
    data.acceptor(idx,:) = data.acceptor(idx,:) - data.total(idx,:)*correctionFactor;
    
    % This version only works for static systems like DNA.
%     stats = traceStat(data);
%     idx = [stats.corr] < -0.5;
%     data.acceptor(idx,:) = data.acceptor(idx,:) - data.total(idx,:)*correctionFactor;
    
    data.recalculateFret();
    data.fileMetadata.ade = correctionFactor;
    
    if ischar(input) && nargout==0
        [p,f,e] = fileparts(input);
        saveTraces( fullfile(p,[f '_ade' e]), data );
    end

end %function


