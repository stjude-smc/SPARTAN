function varargout = separateDays(traceInput,dwtfname)
% separateDays  Split combined traces by date they were obtained.
% 
%    separateDays(FILENAME) splits the traces from a .traces file input into
%    subsets with matching paths in their traceMetadata.ids field. (This field
%    gives the filename of the original .rawtraces file saved in gettraces.)
%    The subsets are saved as "_subsetX.traces" files. If a .dwt file with
%    the same base name is found, "_subsetX.qub.dwt" files are also saved.
%
%    separateDays(TRC_FILE,DWT_FILE) allows a .dwt file to be specified.
%
%    This function is most commonly used to extract traces from distinct days
%    of expeirment from a file that was previously combined from multiple days,
%    for example for calculating error bars.
%    
%    [TRC_SUB,IDL_SUB] = separateDays(...) returns the subsets as a cell array
%    of traces objects and, if applicable, idealization matrices. No output
%    files will be saved.
%
%    See also combineDatasets.

%   Copyright 2016 Cornell University All Rights Reserved.



%% Process input arguments

narginchk(0,2);
nargoutchk(0,3);

% Process traces input
if nargin<1,
    traceInput = getFile;
    if isempty(traceInput), return; end
end
if ischar(traceInput)
    data = loadTraces(traceInput);
    [p,f] = fileparts(traceInput);
    basename = fullfile(p,[f '_']);
elseif isa(traceInput,'TracesFret')
    data = traceInput;
    basename = '';
else
    error('Invalid input');
end

% Look for associated .dwt if not specified.
if nargin<2 && ischar(traceInput),
    dwtfname = findDwt(traceInput);
    dwtfname = dwtfname{1};
end

if ~isempty(dwtfname),
    [dwt,sampling,offsets,model] = loadDWT(dwtfname);
    idl = dwtToIdl(dwt,offsets,data.nFrames,data.nTraces);
end



%%

% Strip out just the path, which we assume is unique to the day of experiments.
ids = {data.traceMetadata.ids};

for i=1:numel(ids),
    idxfile = find( ids{i}=='/'|ids{i}=='\', 1,'last' );
    ids{i} = ids{i}(1:idxfile);
end

% For each unique path, save a file with that subset.
uniquePaths = unique(ids);
nDays = numel(uniquePaths);
[dataSubsets,idlSubsets,idxSubset] = deal( cell(nDays,1) );

for i=1:nDays,
    disp(uniquePaths{i});
    
    % Get subset of traces
    idxSubset{i} = find( strcmp(ids,uniquePaths{i}) );
    dataSubsets{i} = data.getSubset( idxSubset{i} );
    
    if ~isempty(dwtfname),
        idlSubsets{i} = idl( idxSubset{i},:);
    end
    
    % Save files
    if nargout==0,
        outname = sprintf('%ssubset%d.traces',basename,i);
        saveTraces( outname, dataSubsets{i} );
        
        if ~isempty(dwtfname),
            [dwt_subset,offsets_subset] = idlToDwt( idlSubsets{i} );
            outname = strrep(outname,'.traces','.qub.dwt');
            saveDWT( outname, dwt_subset, offsets_subset, model, sampling );
        end
    end
end


output = {dataSubsets,idlSubsets,idxSubset};
[varargout{1:nargout}] = output{1:nargout};




