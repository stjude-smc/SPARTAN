function saveTraces( filename, data )
%SAVETRACES   Save fluorescence/FRET data to a .traces format file.
%
%   SAVETRACES( FNAME, DATA ) saves the data in the Traces object DATA to
%   disk with the file name FNAME (.traces extension). The custom file
%   format mirrors the structure of Traces objects, including metadata.
%
%   See also: loadTraces, TracesFret, TracesFret4, forQuB.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Verify input arguments
narginchk(2,2);
nargoutchk(0,0);

if ~ischar(filename) && ~isa(data,'Traces'),
    error('Invalid input arguments');
end
checkValid(data);

constants = cascadeConstants;
data.fileMetadata(1).software = constants.software;


% Remove existing file (increases speed, for some reason).
if exist(filename,'file'), delete(filename); end


% Write header data
fid = fopen(filename,'Wb');

fwrite( fid, 0,         'uint32' );  %identifies the new traces format.
fwrite( fid, 'TRCS',    'char'   );  %format identifier ("magic")
fwrite( fid, 5,         'uint16' );  %format version number
fwrite( fid, 9,         'uint8'  );  %trace data type (9=single)
fwrite( fid, numel(data.channelNames),    'uint8'  );
fwrite( fid, [data.nTraces data.nFrames], 'uint32' );


% Write channel names (donor, acceptor, fret, etc).
channelNames = strjoin(data.channelNames, char(31));
fwrite( fid, numel(channelNames), 'uint32' );
fwrite( fid, channelNames, 'char' );


% Write fluorescence and FRET data.
fwrite( fid, data.time, 'single' ); %time axis (in seconds)
for i=1:numel(data.channelNames),
    fwrite( fid, data.(data.channelNames{i}), 'single' );
end


% Write metadata.
root.fileMetadata  = data.fileMetadata;
root.traceMetadata = data.traceMetadata;
writeMetadata(fid, root);

fclose(fid);

end  %function saveTraces



function writeMetadata( fid, metadata )
% Writes a metadata object to the currently open file handle fid.

dataTypes = {'char','uint8','uint16','uint32','uint64', ...
                    'int8', 'int16', 'int32', 'int64', ...
             'single','double','logical','cell','struct'};

% Determine data type code and verify it is an allowed type.
type = find( strcmp(class(metadata), dataTypes), 1 );
if isempty(type),
    warning('Skipping invalid metadata type "%s"\n', class(metadata));
    return;
end

szClass = [1  1 2 3 4  1 2 3 4  4 8  1  0 0];  %byte size of primitives
ndim = ndims(metadata);
nBytes = numel(metadata)*szClass(type) +1+4*ndim;


% Write metadata header
fwrite( fid, type-1, 'uint8' );   % field type code, zero-based
fwrite( fid, nBytes, 'uint32');   % content length in bytes
head = ftell(fid);

fwrite( fid, ndim,           'uint8' );   % number of array dimensions (2)
fwrite( fid, size(metadata), 'uint32');   % dimensions of the array


% Write content
if iscell(metadata),
    for i=1:numel(metadata),
        writeMetadata(fid, metadata{i});
    end
    
elseif isstruct(metadata),
    fields = fieldnames(metadata);
    fwrite( fid, numel(fields), 'uint8');   % number of fields
    
    for i=1:numel(fields)
        fname = fields{i};
        fwrite( fid, numel(fname), 'uint8');   % field name length
        fwrite( fid, fname,        'char' );   % field name text
        
        % Try to pack all struct array elements into one array to save
        % space and time. Only works if all elements are the same type.
        isPackable = false;
        firstElement = metadata(1).(fname);
        try
            if ischar(firstElement),
                packed = strjoin({metadata.(fname)}, char(31));
                isPackable = true;
            elseif (isnumeric(firstElement) || islogical(firstElement)),
                packed = cat(ndim+1, metadata.(fname));
                isPackable = size(packed,ndim+1)==numel(metadata) & ...
                             size(packed,ndim+1)>1;
            end
        catch
        end
        fwrite( fid, isPackable, 'uint8' );
        
        % Write the struct array data
        if isPackable,
            writeMetadata( fid, packed );
        else
            for j=1:numel(metadata),
               writeMetadata( fid, metadata(j).(fname) );
            end
        end
    end %for each field in struct
    
else
    fwrite( fid, metadata, class(metadata) );  %primitive types
end


% Record the actual page size of struct/cell types (instead of estimating)
if szClass(type)==0,
    tail = ftell(fid);
    fseek(fid,head-4,-1);
    fwrite( fid, tail-head, 'uint32' );   % actual content length
    fseek(fid,tail,-1);
end


end  %function writeMetadata




