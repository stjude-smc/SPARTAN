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
    outname = fullfile( p, [f '.h5'] );

    % Save channel data to matlab
    data.time = dataIn.time;
    
    for c=1:dataIn.nChannels,
        ch = dataIn.channelNames{c};
        data.(ch) = dataIn.(ch);
    end
    

    dataset = nan(length(dataIn.idxFluor),dataIn.nFrames,dataIn.nTraces);
    dataset(dataIn.idxFluor(1),:,:) = dataIn.donor';
    dataset(dataIn.idxFluor(2),:,:) = dataIn.acceptor';
    
    exposure_frames = dataIn.sampling; % in ms
    
    % Create file 
    fid = H5F.create(outname, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    % H5F_ACC_TRUNC (Flag): Used in H5Fcreate to truncate (overwrite) a file if it already exists, or create it if it does not
    % H5P_DEFAULT (first, file creation): to indicate that the library should use default values for property lists (e.g., file creation or file access property lists
    % H5P_DEFAULT (second, file access): to indicate that the library should use default values for property lists (e.g., file creation or file access property lists

    % Create groups explicitly
    gid_dataset = H5G.create(fid, '/dataset', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    gid_data    = H5G.create(fid, '/dataset/data', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    gid_sources = H5G.create(fid, '/dataset/sources', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    gid_src0    = H5G.create(fid, '/dataset/sources/0', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    
    % Close groups
    H5G.close(gid_src0);
    H5G.close(gid_sources);
    H5G.close(gid_data);
    H5G.close(gid_dataset);
    H5F.close(fid);
    
    % Attributes: /dataset
    h5writeatt(outname, '/dataset', 'format', 'SPARTAN');
    h5writeatt(outname, '/dataset', 'date_created', ...
         string(datetime("today")));
    h5writeatt(outname, '/dataset', 'date_modified', ...
         string(datetime("today")));
    
    h5writeatt(outname, '/dataset/data', 'description', '');
    
    % Dataset: raw
    
    chunkSize = [1 1 1];
    
    h5create(outname, '/dataset/data/raw', size(dataset), ...
        'Datatype', 'double', ...
        'ChunkSize', chunkSize, ...
        'Deflate', 4, ...
        'FillValue', 0);
    
    % Deflate compression level
    % chunksize: non-zero

    h5write(outname, '/dataset/data/raw', dataset);
    
    % Dataset: source_index
    h5create(outname, '/dataset/data/source_index', dataIn.nTraces, ...
        'Datatype', 'int64', ...
        'ChunkSize', dataIn.nTraces, ...
        'Deflate', 4, ...
        'FillValue', int64(0));
    
    % Attributes: /dataset/sources
    h5writeatt(outname, '/dataset/sources', 'source_list', ...
    [p,'/',f,'.traces']);
    
    % Attributes: /dataset/sources/0
    h5writeatt(outname, '/dataset/sources/0', 'source_name', ...
    [p,'/',f,'.traces']);
    display(outname)

 end
 % Create output filename automatically


end %function LoadTracesBinary









