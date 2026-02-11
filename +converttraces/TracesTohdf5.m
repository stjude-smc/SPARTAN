function data = TracesTohdf5( filename, varargin )
% TracesTohd5  Load fluorescence/FRET data from a .traces format  and 
% converts them to hdf5 file.
% For the attributes, I used the attribute type from tmaven
% This file makes .traces compatible with tmaven
% https://github.com/GonzalezBiophysicsLab/tmaven



% Get input filenames from user if not given.
if nargin<1 || isempty(filename),
    filename = getFiles;
end

if ~iscell(filename),
    filename = {filename};
end

data = struct();

for i=1:numel(filename),
    % Load input traces file
    dataIn = loadTraces( filename{i} );
    
    % Create output filename automatically
    [p,f] = fileparts( filename{i} );
    outname = fullfile( p, [f 'tmaven.h5'] );

if exist(outname, 'file') == 2
    try
        delete(outname); % Delete the file
        fprintf('File "%s" deleted successfully.\n', outname);
    catch ME
        % Handle any deletion errors
        fprintf('Error deleting file "%s": %s\n', outname, ME.message);
    end
else
    fprintf('File "%s" does not exist. No action taken.\n', outname);
end

    % Save channel data to matlab
    data.time = dataIn.time;
    
    for c=1:dataIn.nChannels,
        ch = dataIn.channelNames{c};
        data.(ch) = dataIn.(ch);
    end
    

    % Optimize array construction
    nFluor = length(dataIn.idxFluor);
    dataset = zeros(nFluor, dataIn.nFrames, dataIn.nTraces);
    % need [nFluor x nFrames x nTraces], so permute to [nFrames x nTraces] then assign
    dataset(dataIn.idxFluor(1),:,:) = permute(dataIn.donor, [2 1]);
    dataset(dataIn.idxFluor(2),:,:) = permute(dataIn.acceptor, [2 1]);
    
    exposure_frames = dataIn.sampling; % in ms
    
    % Create file strucutre and groups
    fid = H5F.create(outname, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    % H5F_ACC_TRUNC (Flag): Used in H5Fcreate to truncate (overwrite) a file if it already exists, or create it if it does not
    % H5P_DEFAULT (first, file creation): to indicate that the library should use default values for property lists (e.g., file creation or file access property lists
    % H5P_DEFAULT (second, file access): to indicate that the library should use default values for property lists (e.g., file creation or file access property lists

    % Create groups explicitly
    gid_dataset = H5G.create(fid, '/dataset', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    gid_data    = H5G.create(fid, '/dataset/data', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    gid_sources = H5G.create(fid, '/dataset/sources', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    gid_src0    = H5G.create(fid, '/dataset/sources/0', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    
    % Close groups and file (high-level functions will reopen as needed)
    H5G.close(gid_src0);
    H5G.close(gid_sources);
    H5G.close(gid_data);
    H5G.close(gid_dataset);
    H5F.close(fid);
    
    % Batch attribute writes: all attributes written together
    % Attributes: /dataset
    todayStr = string(datetime("today"));
    h5writeatt(outname, '/dataset', 'format', 'SPARTAN');
    h5writeatt(outname, '/dataset', 'date_created', todayStr);
    h5writeatt(outname, '/dataset', 'date_modified', todayStr);
    h5writeatt(outname, '/dataset/data', 'description', '');
    
    % Dataset: raw
    % Chunk size is improtant for efficient writes
    % Typical good chunk sizes: 100-1000 elements per dimension
    nFluor = size(dataset, 1);
    nFrames = size(dataset, 2);
    nTraces = size(dataset, 3);
    chunkSize = [min(100, nFluor), min(1000, nFrames), min(1000, nTraces)];
    
    h5create(outname, '/dataset/data/raw', size(dataset), ...
        'Datatype', 'double', ...
        'ChunkSize', chunkSize, ...
        'Deflate', 2, ...  % compression level 
        'FillValue', 0);
    
    % Deflate: compression level
    % chunksize: non-zero values for saving data in chunks

    h5write(outname, '/dataset/data/raw', dataset);
    
    % Dataset: source_index
    % Optimized: use reasonable chunk size (entire array if small, otherwise cap it)
    sourceIndexChunkSize = min(10000, max(100, dataIn.nTraces));
    h5create(outname, '/dataset/data/source_index', dataIn.nTraces, ...
        'Datatype', 'int64', ...
        'ChunkSize', sourceIndexChunkSize, ...
        'Deflate', 2, ...  % compression level 
        'FillValue', int64(0));
    
    % Attributes: /dataset/sources
    h5writeatt(outname, '/dataset/sources', 'source_list', ...
    [p,'/',f,'.traces']);
    
    % Attributes: /dataset/sources/0
    h5writeatt(outname, '/dataset/sources/0', 'source_name', ...
    [p,'/',f,'.traces']);
    display(outname)

 end


end %function LoadTracesBinary









