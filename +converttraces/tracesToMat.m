function tracesToMat( filenames )
% tracesToMat  Convert .traces files to .mat format (native Matlab)

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% 
% % Get input filenames from user if not given.
% if nargin<1 || isempty(filenames),
%     filenames = getFiles;
% end
% 
if ~iscell(filenames),
    filenames = {filenames};
end


for i=1:numel(filenames),

    dataTypes = {'char','uint8','uint16','uint32','uint64', ...
                    'int8', 'int16', 'int32', 'int64', ...
             'single','double','logical','cell','struct'};  %zero-based

    % Create output filename automatically
    [p,f] = fileparts( filenames{i} );
    outname = fullfile( p, [f '.mat'] );

    fid = fopen(filenames{i},'r');

    fread(fid,1,'uint32');        % ignore version
    fread(fid,[1,4],'*char');     % ignore signature
    fread(fid,1,'uint16');        % ignore format version

    dataType  = fread( fid, 1, '*uint8' );   % data class (9=single)
    nChannels = fread( fid, 1, '*uint8' );
    nTraces   = fread( fid, 1, 'uint32' );
    traceLen  = fread( fid, 1, 'uint32' );

    szNames      = fread( fid, 1, 'uint32' );
    channelNames = fread( fid, [1,szNames], '*char' );
    channelNames = strsplit(channelNames,char(31), 'CollapseDelimiters',false);

    type = dataTypes{dataType+1};
    data.time = fread( fid, [1,traceLen], type );  %time axis (s)
    
    indexes = 1:nTraces;
    useAll = isequal(indexes, 1:nTraces);

    if dataType==9,
    star = '*';
    else
        disp('there is a problem')
    end

    for i=1:nChannels,
        d = fread( fid, [nTraces,traceLen], [star type] );
        
        if useAll,  %faster if indexing is avoided.
            data.(channelNames{i}) = d;
        else
            data.(channelNames{i}) = d(indexes,:);
        end
    end

    fclose(fid);
    
    % Save to file.
    save(outname,'data');
end


end