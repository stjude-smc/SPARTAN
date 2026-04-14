function tracesToNpz( filenames,deeplasidatapath,pythonExepath )
% tracesToNpz  Convert .traces files to .npz format (native Python numpy array)

%  I used the logic in tracesToMat

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
    outname = fullfile( deeplasidatapath, [f '.mat'] );
    outname1 = fullfile( deeplasidatapath, [f 'dlasi.npz'] );

    fid = fopen(filenames{i},'r');

    fread(fid,1,'uint32');        % ignore version
    fread(fid,[1,4],'*char');     % ignore signature
    fread(fid,1,'uint16');        % ignore format version

    dataType  = fread( fid, 1, '*uint8' );   % data class (9=single)
    nChannels = fread( fid, 1, '*uint8' )-3; %'donor','acceptor','acceptor2','acceptorDirect','acceptorDirect2','stoichiometry','fret','fret2'
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

    arr_0 = nan(nTraces,traceLen,nChannels+1); % 1)donor, 2)acceptor,3)accpetor2,4)acceptorDirect,5)accptorDirect2 and 6)acceptor2Direct(not existing)
    for i=1:nChannels,
        d = fread( fid, [nTraces,traceLen], [star type] );
        
        if useAll,  %faster if indexing is avoided.
            data.(channelNames{i}) = d;
            arr_0(:,:,i) = d;
        else
            data.(channelNames{i}) = d(indexes,:);
            arr_0(:,:,i) =  d(indexes,:);
        end
    end

    arr_0(:,:,6) =zeros(size(d));

    fclose(fid);
    
    % Save to file.
    save(outname,'arr_0');
    

thisDir = fileparts(mfilename('fullpath'));   % .../SPARTAN/+converttraces
pythonExe = pythonExepath; %'/xxx/venv/bin/python';  % user or GUI
matPath = outname ; %'/xxx/deeplasi-main/functions/deeplearning/data/aarondata.mat';
npzPath = outname1;%'/xxx/deeplasi-main/functions/deeplearning/data/aarondata_npz.npz';

if ispc
    script    =  fullfile(thisDir, 'tracesMat2Npz.py');
    cmd = sprintf('"%s" "%s" "%s" "%s"', ...
                pythonExe, ...
                script, ...
                matPath, ...
                npzPath);
elseif isunix 
    script    =  fullfile(thisDir, 'run_python.sh');
    cmd = sprintf('bash "%s" "%s" "%s" "%s"',script, pythonExe,matPath,npzPath);
end



[status, out] = system(cmd);
disp(out);

if status ~= 0
    error('MAT to NPZ conversion failed (status=%d):\n%s', status, out);
end


end


end

