function applyLocalGamma(ffcorr,files)
% applyLocalGamma  correction for local variance in detection efficiency.
%
%   applyLocalGamma(FILES) where FILES the path (char array) to a
%   .traces file, or a cell array of scuh paths, calculates a spatial
%   mapping function R(x,y) of the acceptor channel detection sensitivity
%   relative to the donor channel.
%   FILES should be obtained from movies with a standard sample,
%   such as DNA oligos, with a fixed FRET efficiency of ~0.5.
%   R(x,y) can be used by applyLocalGamma for flat-field correction to
%   correct for non-uniform detection efficiency.
%   
%   applyLocalGamma() without input arguments will prompt for files.
%   applyLocalGamma() without input arguments will prompt for files.
% 
% See also: calcLocalGamma, calc_gamma, gammacorrect, correctTraces.

% Copyright 2025 St Jude Children's Research Hospital. All Rights Reserved.

% FIXME:
%  - does not adjust gamma value in traceMetadata -- it should!
%  - correction of cropped (high time resolution) movies is untested.
%  - consider returning corrected data directly if input is Traces object.


% Select input files
narginchk(0,2);

if nargin<1 || isempty(ffcorr)
    ffcorr = getFile('*.mat','Select correction file');
    if isempty(ffcorr), return; end  %user hit cancel.
end
if ischar(ffcorr) || isstring(ffcorr)
    ffcorr = load(ffcorr);
end
assert( isstruct(ffcorr) && all(isfield(ffcorr,{'result','nX','nY'})), ...
        'Invalid flat-field correction calibration file.' );

if nargin<2 || isempty(files)
    filter = {'*.traces','Binary Traces Files (*.traces)'; ...
              '*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
              '*.*','All Files (*.*)'};
    files = getFiles(filter,'Select traces to correct',false);
    if isempty(files), return; end  %user hit cancel
end


% Apply flat field correction to each files
for i=1:numel(files)
    data = loadTraces(files{i});

    if ~all(isfield(data.fileMetadata,{'nX','nY'}))
        error('SPARTAN:applyLocalGamma:fileMetadataMissing', ...
                'Traces input missing field size; assuming same as correction file.');
    end
    assert(data.fileMetadata.nX==ffcorr.nX, 'Field of view mismatch');
    
    x = [data.traceMetadata.donor_x]';
    y = [data.traceMetadata.donor_y]' + ffcorr.nY-data.fileMetadata.nY;
    data.acceptor = data.acceptor ./ ffcorr.result(x,y);
    data.recalculateFret();

    % Adjust gamma metadata to match -- FIXME, untested
    idxAcc = find(strcmpi(data.channelNames,'acceptor'),1);
    for t=1:data.nTraces
       data.traceMetadata(t).scaleFluor(idxAcc) = data.traceMetadata(t).scaleFluor(idxAcc) ...
           / ffcorr.result(x(t),y(t));
    end

    [p,f,e] = fileparts(files{i});
    data.save( fullfile(p,[f '_ffcorr' e]) );
end


end %function


