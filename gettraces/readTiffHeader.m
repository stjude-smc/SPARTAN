function [stack,dataOffsets] = readTiffHeader(filename,indices)
% readTiffHeader    Parses a TIFF file (directory structure)
% This function is designed for use with the Movie_TIFF class. It should load
% all TIFF-based formats, including MetaMorph STK, MetaMorph TIFF, Zeiss LSM,
% and traditional TIFF stacks.
%
%    
%
% WARNING: these scripts have only been tested on MetaMorph STK files
% (old format) with 16-bit greyscale data and may not work at all with other
% formats.
%
% This script is based on the "tiffread29" code (May 10, 2010) by Francois
% Nedelec (EMBL) under the GNU General Public License (version 3). This code
% is also available under the same license (see the LICENSE file).
% See http://www.embl.de/~nedelec/misc/index.html
% 
% The original code was modified significantly to improve performance at loading
% MetaMorph STK files:
% 1) Assumes litte-endian bit order. Not compatible with PowerPC, Solaris, etc.
% 2) Removed wait-bar.
% 3) Optimize reading of MM-specific info and not storing the raw text in every
%    plane (as was the case before).
% 4) Assumes strips always completely consolidate into one plane!
%      -- DT (02/16/08)
% QUESTION: where is the stack free-hand annotation??

%   Copyright 2007-2015 Cornell University All Rights Reserved.


%Optimization: join adjacent TIF strips: this results in faster reads
consolidateStrips = 1;

%without argument, we ask the user to choose a file:
if nargin < 1
    [filename, pathname] = uigetfile('*.tif;*.stk;*.lsm', 'select image file');
    filename = [ pathname, filename ];
end

if (nargin<=1);  indices = []; end


% not all valid tiff tags have been included, as they are really a lot...
% if needed, tags can easily be added to this code
% See the official list of tags:
% http://partners.adobe.com/asn/developer/pdfs/tn/TIFF6.pdf
%
% the structure IMG is returned to the user, while TIF is not.
% so tags usefull to the user should be stored as fields in IMG, while
% those used only internally can be stored in TIF.

global TIF;
TIF = [];

%counters for the number of images read and skipped
img_skip  = 0;
img_read  = 1;

%% set defaults values :
TIF.SampleFormat     = 1;
TIF.SamplesPerPixel  = 1;
TIF.BOS              = 'ieee-le';          %byte order string

if  isempty(findstr(filename,'.'))
    filename = [filename,'.tif'];
end

TIF.file = fopen(filename,'r','l');
if TIF.file == -1
    stkname = strrep(filename, '.tif', '.stk');
    TIF.file = fopen(stkname,'r','l');
    if TIF.file == -1
        error(['File "',filename,'" not found.']);
    else
        filename = stkname;
    end
end
[s, m] = fileattrib(filename);

% obtain the full file path:
filename = m.Name;

% find the file size in bytes:
% m = dir(filename);
% filesize = m.bytes;


%% read header
% read byte order: II = little endian, MM = big endian
byte_order = fread(TIF.file, 2, '*char');
if ( strcmp(byte_order', 'II') )
    TIF.BOS = 'ieee-le';                                % Intel little-endian format
elseif ( strcmp(byte_order','MM') )
    TIF.BOS = 'ieee-be';
    error('Big-endian byte-ordering not supported.');
else
    error('This is not a TIFF file (no MM or II).');
end


%% ---- read in a number which identifies file as TIFF format
tiff_id = fread(TIF.file,1,'uint16');
if (tiff_id ~= 42)
    error('This is not a TIFF file (missing 42).');
end

%% ---- read the byte offset for the first image file directory (IFD)
TIF.img_pos = fread(TIF.file, 1, 'uint32');

while  TIF.img_pos ~= 0 

    clear IMG;
    IMG.filename = filename;
    % move in the file to the first IFD
    status = fseek(TIF.file, TIF.img_pos, -1);
    if status == -1
        error('invalid file offset (error on fseek)');
    end

    %disp(strcat('reading img at pos :',num2str(TIF.img_pos)));

    %read in the number of IFD entries
    num_entries = fread(TIF.file,1,'uint16');
    %disp(strcat('num_entries =', num2str(num_entries)));

    %read and process each IFD entry
    for i = 1:num_entries

        % save the current position in the file
        file_pos  = ftell(TIF.file);

        % read entry tag
        TIF.entry_tag = fread(TIF.file, 1, 'uint16');
        % read entry
        entry = readIFDentry;
        %disp(strcat('reading entry <',num2str(TIF.entry_tag),'>'));

        switch TIF.entry_tag
            case 254
                TIF.NewSubfiletype = entry.val;
            case 256         % image width - number of column
                IMG.width          = entry.val;
            case 257         % image height - number of row
                IMG.height         = entry.val;
                TIF.ImageLength    = entry.val;
            case 258         % BitsPerSample per sample
                TIF.BitsPerSample  = entry.val;
                TIF.BytesPerSample = TIF.BitsPerSample / 8;
                IMG.bits           = TIF.BitsPerSample(1);
                %fprintf('BitsPerSample %i %i %i\n', entry.val);
            case 259         % compression
                if ( entry.val ~= 1 )
                    error(['Compression format ', num2str(entry.val),' not supported.']);
                end
            case 262         % photometric interpretation
                TIF.PhotometricInterpretation = entry.val;
                if ( TIF.PhotometricInterpretation == 3 )
                    warning('tiffread2:LookUp', 'Ignoring TIFF look-up table');
                end
            case 266    % FillOrder
                assert( entry.val==1, 'FillOrder=2 not supported' );
            case 269
                IMG.document_name  = entry.val;
            case 270         % comments:
                IMG.info           = entry.val;
            case 271
                IMG.make           = entry.val;
            case 273         % strip offset
                TIF.StripOffsets   = entry.val;
                TIF.StripNumber    = entry.cnt;
                %fprintf('StripNumber = %i, size(StripOffsets) = %i %i\n', TIF.StripNumber, size(TIF.StripOffsets));
            case 277         % sample_per pixel
                TIF.SamplesPerPixel  = entry.val;
                %fprintf('Color image: sample_per_pixel=%i\n',  TIF.SamplesPerPixel);
            case 278         % rows per strip
                TIF.RowsPerStrip   = entry.val;
            case 279         % strip byte counts - number of bytes in each strip after any compressio
                TIF.StripByteCounts= entry.val;
            case 282         % X resolution
                IMG.x_resolution   = entry.val;
            case 283         % Y resolution
                IMG.y_resolution   = entry.val;
            case 284         %planar configuration describe the order of RGB
                TIF.PlanarConfiguration = entry.val;
            case 296         % resolution unit
                IMG.resolution_unit= entry.val;
            case 305         % software
                IMG.software       = entry.val;
            case 306         % datetime
                IMG.datetime       = entry.val;
            case 315
                IMG.artist         = entry.val;
            case 317        %predictor for compression
                if (entry.val ~= 1); error('unsuported predictor value'); end
            case 320         % color map
                IMG.cmap           = entry.val;
                IMG.colors         = entry.cnt/3;
            case 339
                TIF.SampleFormat   = entry.val;
            case 33550       % GeoTIFF ModelPixelScaleTag
                IMG.ModelPixelScaleTag    = entry.val;
            case 33628       %metamorph specific data
                IMG.MM_private1    = entry.val;
            case 33629       %this tag identify the image as a Metamorph stack!
                TIF.MM_stack       = entry.val;
                TIF.MM_stackCnt    = entry.cnt;
            case 33630       %metamorph stack data: wavelength
                TIF.MM_wavelength  = entry.val;
            case 33631       %metamorph stack data: gain/background?
                TIF.MM_private2    = entry.val;
            case 33922       % GeoTIFF ModelTiePointTag
                IMG.ModelTiePointTag    = entry.val;
            case 34412       % Zeiss LSM data
                LSM_info           = entry.val;
            case 34735       % GeoTIFF GeoKeyDirectory
                IMG.GeoKeyDirTag       = entry.val;
            case 34737       % GeoTIFF GeoASCIIParameters
                IMG.GeoASCII       = entry.val;
            case 42113       % GeoTIFF GDAL_NODATA
                IMG.GDAL_NODATA    = entry.val;
            otherwise
                fprintf( 'Ignored TIFF entry with tag %i (cnt %i)\n', TIF.entry_tag, entry.cnt);
        end
        
        % calculate bounding box  if you've got the stuff
        %if isfield(IMG, 'ModelPixelScaleTag') && isfield(IMG, 'ModelTiePointTag') && isfield(IMG, 'height')&& isfield(IMG, 'width'),
        %    IMG.North=IMG.ModelTiePointTag(5)-IMG.ModelPixelScaleTag(2)*IMG.ModelTiePointTag(2);
        %    IMG.South=IMG.North-IMG.height*IMG.ModelPixelScaleTag(2);
        %    IMG.West=IMG.ModelTiePointTag(4)+IMG.ModelPixelScaleTag(1)*IMG.ModelTiePointTag(1);
        %    IMG.East=IMG.West+IMG.width*IMG.ModelPixelScaleTag(1);
        %end

        % move to next IFD entry in the file
        status = fseek(TIF.file, file_pos+12, -1);
        if status == -1
            error('invalid file offset (error on fseek)');
        end
    end

    %Planar configuration is not fully supported
    %Per tiff spec 6.0 PlanarConfiguration irrelevent if SamplesPerPixel==1
    %Contributed by Stephen Lang
    if (TIF.SamplesPerPixel ~= 1) && ( ~isfield(TIF, 'PlanarConfiguration') || TIF.PlanarConfiguration == 1 )
        error('PlanarConfiguration = 1 is not supported');
    end

    %total number of bytes per image:
    PlaneBytesCnt = IMG.width * IMG.height * TIF.BytesPerSample;

    %% try to consolidate the TIFF strips if possible
    
    if consolidateStrips
        %Try to consolidate the strips into a single one to speed-up reading:
        BytesCnt = TIF.StripByteCounts(1);

        if BytesCnt < PlaneBytesCnt

            ConsolidateCnt = 1;
            %Count how many Strip are needed to produce a plane
            while TIF.StripOffsets(1) + BytesCnt == TIF.StripOffsets(ConsolidateCnt+1)
                ConsolidateCnt = ConsolidateCnt + 1;
                BytesCnt = BytesCnt + TIF.StripByteCounts(ConsolidateCnt);
                if ( BytesCnt >= PlaneBytesCnt ); break; end
            end

            %Consolidate the Strips
            if ( BytesCnt <= PlaneBytesCnt(1) ) && ( ConsolidateCnt > 1 )
                %fprintf('Consolidating %i stripes out of %i', ConsolidateCnt, TIF.StripNumber);
                TIF.StripByteCounts = [BytesCnt; TIF.StripByteCounts(ConsolidateCnt+1:TIF.StripNumber ) ];
                TIF.StripOffsets = TIF.StripOffsets( [1 , ConsolidateCnt+1:TIF.StripNumber] );
                TIF.StripNumber  = 1 + TIF.StripNumber - ConsolidateCnt;
            end
        end
        
        assert( numel(TIF.StripOffsets)==1, 'Failed to conslidate strips' );
    end

    %% read the next IFD address:
    TIF.img_pos = fread(TIF.file, 1, 'uint32');
    %if (TIF.img_pos) disp(['next ifd at', num2str(TIF.img_pos)]); end

    if isfield( TIF, 'MM_stack' )

        if isempty(indices),
            indices = 1:TIF.MM_stackCnt;
        else
            indices = indices( indices<=TIF.MM_stackCnt );
        end
        
        % Save MetaMorph metadata
        if isfield(IMG,'info'),
            %MM = parseMetamorphInfo(IMG.info, TIF.MM_stackCnt);
            MM = parseMetamorphInfo(IMG.info, 1);
        else
            MM = [];
        end
        
        IMG = rmfield( IMG, 'info' );

        %this loop reads metamorph stacks:
        % NOTE: only some of the image info is really available on a per-frame
        % basis: MM_stack, MM_wavelength, MM_private2, MM. that's it. There's
        % some major wasted space duplicating the others.
        imgByteSize = IMG.width*IMG.height*TIF.BytesPerSample;
        %dataOffsets = zeros(1,numel(indices));
        
        dataOffsets = TIF.StripOffsets+(indices-1).*imgByteSize;
        
        ii=1;
%         for ii = indices

            TIF.StripCnt = 1;
            
            % Add MetaMorph metadata to the image.
            [ IMG.MM_stack, IMG.MM_wavelength, IMG.MM_private2 ] = splitMetamorph(ii);
            if ~isempty(MM),
                IMG.MM = MM(img_read);
            end
            
            % Save data location.
            % WARNING: assumes consolidated strips and single-channel data.
            dataOffsets(img_read) = TIF.StripOffsets+(ii-1)*imgByteSize;
            
            % Save the current image into the stack.
            stack(img_read) = IMG;
            %img_read = img_read + 1;
%         end
        
        
        break;

    else

        %this part reads a normal TIFF stack:
        
        read_img = any( img_skip+img_read == indices );
        if exist('stack','var')
            if IMG.width ~= stack(1).width || IMG.height ~= stack(1).height
                %setting read_it=0 will skip dissimilar images:
                %comment-out the line below to allow dissimilar stacks
                read_img = 0;
            end
        end
        
        if read_img
            TIF.StripCnt = 1;
            %read the image channels
            %for c = 1:TIF.SamplesPerPixel
            %    IMG.data{c} = read_plane(0, IMG.width, IMG.height, c);
            %end
            dataOffsets(img_read) = TIF.StripOffsets;

            try
                stack(img_read) = IMG;  % = orderfields(IMG);
                img_read = img_read + 1;
            catch
                fprintf('Tiffread skipped dissimilar image %i\n', img_read+img_skip);
                img_skip = img_skip + 1;
             end
            
            if  all( img_skip+img_read > indices )
                break;
            end

        else
            img_skip = img_skip + 1;
        end

    end
end


%% duplicate the LSM info
if exist('LSM_info', 'var')
    for i = 1:numel(stack)
        stack(i).lsm = LSM_info;
    end
end


%% return

if ~ exist('stack', 'var')
    stack = struct([]);
end

%clean-up
fclose(TIF.file);



end




%% ==================sub-functions that reads an IFD entry:===================


function [nbBytes, matlabType] = convertType(tiffType)
switch (tiffType)
    case 1
        nbBytes=1;
        matlabType='uint8';
    case 2
        nbBytes=1;
        matlabType='uchar';
    case 3
        nbBytes=2;
        matlabType='uint16';
    case 4
        nbBytes=4;
        matlabType='uint32';
    case 5
        nbBytes=8;
        matlabType='uint32';
    case 7
        nbBytes=1;
        matlabType='uchar';
    case 11
        nbBytes=4;
        matlabType='float32';
    case 12
        nbBytes=8;
        matlabType='float64';
    otherwise
        error('tiff type %i not supported', tiffType)
end
end

%% ==================sub-functions that reads an IFD entry:===================

function  entry = readIFDentry()

global TIF;
entry.tiffType = fread(TIF.file, 1, 'uint16');
entry.cnt      = fread(TIF.file, 1, 'uint32');
%disp(['tiffType =', num2str(entry.tiffType),', cnt = ',num2str(entry.cnt)]);

[ entry.nbBytes, entry.matlabType ] = convertType(entry.tiffType);

if entry.nbBytes * entry.cnt > 4
    %next field contains an offset:
    offset = fread(TIF.file, 1, 'uint32');
    %disp(strcat('offset = ', num2str(offset)));
    status = fseek(TIF.file, offset, -1);
    if status == -1
        error('invalid file offset (error on fseek)');
    end

end


if TIF.entry_tag == 33629   % metamorph 'rationals'
    entry.val = fread(TIF.file, 6*entry.cnt, entry.matlabType);
elseif TIF.entry_tag == 34412  %TIF_CZ_LSMINFO
    entry.val = readLSMinfo;
else
    if entry.tiffType == 5
        entry.val = fread(TIF.file, 2*entry.cnt, entry.matlabType);
    else
        entry.val = fread(TIF.file, entry.cnt, entry.matlabType);
    end
end

if ( entry.tiffType == 2 );
    entry.val = char(entry.val');
end

end


%% =============distribute the metamorph infos to each frame:
function [MMstack, MMwavelength, MMprivate2] = splitMetamorph(imgCnt)

global TIF;

MMstack = [];
MMwavelength = [];
MMprivate2 = [];

if TIF.MM_stackCnt == 1
    return;
end

left  = imgCnt - 1;

if isfield( TIF, 'MM_stack' )
    S = length(TIF.MM_stack) / TIF.MM_stackCnt;
    MMstack = TIF.MM_stack(S*left+1:S*left+S);
end

if isfield( TIF, 'MM_wavelength' )
    S = length(TIF.MM_wavelength) / TIF.MM_stackCnt;
    MMwavelength = TIF.MM_wavelength(S*left+1:S*left+S);
end

if isfield( TIF, 'MM_private2' )
    S = length(TIF.MM_private2) / TIF.MM_stackCnt;
    MMprivate2 = TIF.MM_private2(S*left+1:S*left+S);
end

end


%% %%  Parse the Metamorph camera info tag into respective fields
% EVBR 2/7/2005, FJN Dec. 2007
function mm = parseMetamorphInfo(info, cnt)

% info   = regexprep(info, '\r\n|\o0', '\n');
info( info==char(0)  ) =char(10);
info( info==char(13) ) =char(10);

parse  = textscan(info, '%s %s', 'Delimiter', ':');
tokens = parse{1};
values = parse{2};

first = char(tokens(1,1));

k = 0;
mm = struct('Exposure', zeros(cnt,1));
for i=1:size(tokens,1)
    tok = char(tokens(i,1));
    val = char(values(i,1));
    %fprintf( '"%s" : "%s"\n', tok, val);
    if strcmp(tok, first)
        k = k + 1;
        if k>cnt, return; end
    end
    if strcmp(tok, 'Exposure')
        [v, c, e, pos] = sscanf(val, '%i');
        unit = val(pos:length(val));
        %return the exposure in milli-seconds
        switch( unit )
            case 'ms'
                mm(k).Exposure = v;
            case 'msec'
                mm(k).Exposure = v;
            case 's'
                mm(k).Exposure = v * 1000;
            otherwise
                warning('tiffread2:Unit', ['Exposure unit "',unit,'" not recognized']);
                mm(k).Exposure = v;
        end
    else
        switch tok
            case 'Binning'
                % Binning: 1 x 1 -> [1 1]
                mm(k).Binning = sscanf(val, '%d %*s %d')';
            case 'Region'
                mm(k).Region = sscanf(val, '%d %*s %d, offset at (%d, %d)')';
            otherwise
                field = strrep(tok,' ','_');
                
                if strcmp(val, 'Off')
                    mm(k).(field)=0;
                elseif strcmp(val, 'On')
                    mm(k).(field)=1;
                elseif isstrprop(val,'digit')
                    mm(k).(field)=str2double(val);
                else
                    mm(k).(field)=val;
                end
        end
    end
end

end

%% ==============partial-parse of LSM info:

function R = readLSMinfo()

% Read part of the LSM info table version 2
% this provides only very partial information, since the offset indicate that
% additional data is stored in the file
global TIF;

R.MagicNumber            = sprintf('0x%09X',fread(TIF.file, 1, 'uint32'));
StructureSize          = fread(TIF.file, 1, 'uint32');
R.DimensionX             = fread(TIF.file, 1, 'uint32');
R.DimensionY             = fread(TIF.file, 1, 'uint32');
R.DimensionZ             = fread(TIF.file, 1, 'uint32');
R.DimensionChannels      = fread(TIF.file, 1, 'uint32');
R.DimensionTime          = fread(TIF.file, 1, 'uint32');
R.IntensityDataType      = fread(TIF.file, 1, 'uint32');
R.ThumbnailX             = fread(TIF.file, 1, 'uint32');
R.ThumbnailY             = fread(TIF.file, 1, 'uint32');
R.VoxelSizeX             = fread(TIF.file, 1, 'float64');
R.VoxelSizeY             = fread(TIF.file, 1, 'float64');
R.VoxelSizeZ             = fread(TIF.file, 1, 'float64');
R.OriginX                = fread(TIF.file, 1, 'float64');
R.OriginY                = fread(TIF.file, 1, 'float64');
R.OriginZ                = fread(TIF.file, 1, 'float64');
R.ScanType               = fread(TIF.file, 1, 'uint16');
R.SpectralScan           = fread(TIF.file, 1, 'uint16');
R.DataType               = fread(TIF.file, 1, 'uint32');
OffsetVectorOverlay    = fread(TIF.file, 1, 'uint32');
OffsetInputLut         = fread(TIF.file, 1, 'uint32');
OffsetOutputLut        = fread(TIF.file, 1, 'uint32');
OffsetChannelColors    = fread(TIF.file, 1, 'uint32');
R.TimeInterval           = fread(TIF.file, 1, 'float64');
OffsetChannelDataTypes = fread(TIF.file, 1, 'uint32');
OffsetScanInformation  = fread(TIF.file, 1, 'uint32');
OffsetKsData           = fread(TIF.file, 1, 'uint32');
OffsetTimeStamps       = fread(TIF.file, 1, 'uint32');
OffsetEventList        = fread(TIF.file, 1, 'uint32');
OffsetRoi              = fread(TIF.file, 1, 'uint32');
OffsetBleachRoi        = fread(TIF.file, 1, 'uint32');
OffsetNextRecording    = fread(TIF.file, 1, 'uint32');

% There are more information stored in this table, which is not read here


%read real acquisition times:
if ( OffsetTimeStamps > 0 )
    
    status = fseek(TIF.file, OffsetTimeStamps, -1);
    if status == -1
        error('error on fseek');
    end
    
    StructureSize          = fread(TIF.file, 1, 'int32');
    NumberTimeStamps       = fread(TIF.file, 1, 'int32');
    for i=1:NumberTimeStamps
        R.TimeStamp(i)       = fread(TIF.file, 1, 'float64');
    end
    
    %calculate elapsed time from first acquisition:
    R.TimeOffset = R.TimeStamp - R.TimeStamp(1);
    
end


end

