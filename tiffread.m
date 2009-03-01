function [stack, time] = tiffreadHDK(filename, img_first, img_last)
% tiffread, version 2-HDK
% WARNING: this is a customized version and has only been tested on 16-bit
% greyscale stacks generated by MetaMorph.  Most other functionality has
% been left intact, but we do not know if it still works.
%
% [stack, nbImages] = tiffread;
% [stack, nbImages] = tiffread(filename);
% [stack, nbImages] = tiffread(filename, imageIndex);
% [stack, nbImages] = tiffread(filename, firstImageIndex, lastImageIndex);
%
% Reads 8,16,32 bits uncompressed grayscale and (some) color tiff files,
% as well as stacks or multiple tiff images, for example those produced
% by metamorph or NIH-image. However, the entire TIFF standard is not
% supported (but you may extend it).
%
% The function can be called with a file name in the current directory,
% or without argument, in which case it pop up a file openning dialog
% to allow manual selection of the file.
% If the stacks contains multiples images, loading can be restricted by
% specifying the first and last images to read, or just one image to read.
%
% at return, nbimages contains the number of images read, and S is a vector
% containing the different images with some additional informations. The
% image pixels values are stored in the field .data, for gray level images,
% or in the fields .red, .green and .blue
% the pixels values are in the native (integer) format,
% and must be converted to be used in most matlab functions.
%
% Example:
% im = tiffread('spindle.stk');
% imshow( double(im(5).data), [] );
%
% Francois Nedelec, EMBL, Copyright 1999-2004.
% Last modified July 7th, 2004 at Woods Hole during the physiology course.
% Thanks to Kendra Burbank for suggesting the waitbar
%
% Please, help us improve this software: send us feedback/bugs/suggestions
% This software is provided at no cost by a public research institution.
% However, postcards are always welcome!
%
% Francois Nedelec
% nedelec (at) embl.de
% Cell Biology and Biophysics, EMBL; Meyerhofstrasse 1; 69117 Heidelberg; Germany
% http://www.embl.org
% http://www.cytosim.org
%
% Slight changes made to read *.lsm files produced by Zeiss LSM 510 Confocal Laser Scanning Microscopes
% by Jan-Ulrich Kreft (kreft@uni-bonn.de), 08/06/2005
%
% Assume little-edian bit ordering for host and data files = major
% performance increase; other performance improvements; removed waitbar (no
% longer needed).  Removed color TIFF support for simplicity.
% HACK: assumes strips always completely consolidate into one plane!!!
% NOTE: consolidating plane reads into one giant fread provides no
% perfomance benefit because dimensions of the resulting 3D matrix
% must be permuted (file is row major, Matlab is column major).
% 02/16/08

%Optimization: join adjacent TIF strips: this results in faster reads
consolidateStrips = true;

%if there is no argument, we ask the user to choose a file:
if (nargin == 0)
    [filename, pathname] = uigetfile('*.tif;*.stk;*.lsm', 'select image file');
    filename = [ pathname, filename ];
end

if (nargin<=1)
    img_first = 1; 
    img_last = 10000; 
end
if (nargin==2)  
    img_last = img_first;            
end


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
img_read  = 0;

% set defaults values :
TIF.SampleFormat     = 1;
TIF.SamplesPerPixel  = 1;
TIF.maxres           = 1; % Jan: highest resolution of images in lsm file

if  isempty(findstr(filename,'.'))
    filename = [filename,'.tif'];
end

[TIF.file, message] = fopen(filename,'r');
if TIF.file == -1
    filename = strrep(filename, '.tif', '.stk');
    [TIF.file, message] = fopen(filename,'r');
    if TIF.file == -1
        error(['file <',filename,'> not found.']);
    end
end


% read header
% read byte order: II = little endian, MM = big endian
%byte_order = setstr(fread(TIF.file, 2, 'uchar'));

byte_order =fread(TIF.file, 2, '*char'); %NOC: 9-23-05

if ( strcmp(byte_order','II') )
    %Intel PC format - ok
elseif ( strcmp(byte_order','MM') )
    error('Big-endian byte-ordered files not supported');
else
    error('This is not a TIFF file (no MM or II).');
end

%----- read in a number which identifies file as TIFF format
tiff_id = fread(TIF.file,1,'uint16');
if (tiff_id ~= 42)
    error('This is not a TIFF file (missing 42).');
end

%----- read the byte offset for the first image file directory (IFD)
ifd_pos = fread(TIF.file,1,'uint32');

while (ifd_pos ~= 0)  %follow directory pointers until end

    clear IMG;
    IMG.filename = [pwd filesep filename];
    % move in the file to the first IFD
    fseek(TIF.file, ifd_pos, -1);
    %disp(strcat('reading img at pos :',num2str(ifd_pos)));

    %read in the number of IFD entries
    num_entries = fread(TIF.file,1,'uint16');
    %disp(strcat('num_entries =', num2str(num_entries)));

    %read and process each IFD entry
    for i = 1:num_entries

        % save the current position in the file
        file_pos  = ftell(TIF.file);

        % read entry tag
        TIF.entry_tag = fread(TIF.file, 1, 'uint16');
        entry = readIFDentry;
        %disp(strcat('reading entry <',num2str(TIF.entry_tag),'>'));

        switch TIF.entry_tag
            case 254
                TIF.NewSubfiletype = entry.val;
            case 256         % image width - number of column
                IMG.width          = entry.val;
                Rwidth = IMG.width;
                if entry.val > TIF.maxres; TIF.maxres = entry.val; end % Jan: set maxres to the highest resolution occuring in the lsm file
            case 257         % image height - number of row
                IMG.height         = entry.val;
                Rheight             = IMG.height;
                TIF.ImageLength    = entry.val;
            case 258         % BitsPerSample per sample
                TIF.BitsPerSample  = entry.val;
                TIF.BytesPerSample = TIF.BitsPerSample / 8;
                IMG.bits           = TIF.BitsPerSample(1);
                %fprintf(1,'BitsPerSample %i %i %i\n', entry.val);
            case 259         % compression
                if (entry.val ~= 1), error('Compression format not supported.'); end;
            case 262         % photometric interpretation
                TIF.PhotometricInterpretation = entry.val;
                if ( TIF.PhotometricInterpretation == 3 )
                    fprintf(1, 'warning: ignoring the look-up table defined in the TIFF file');
                end
            case 269
                IMG.document_name  = entry.val;
            case 270         % comment:
                TIF.info           = entry.val;
            case 271
                IMG.make           = entry.val;
            case 273         % strip offset
                TIF.StripOffsets   = entry.val;
                TIF.StripNumber    = entry.cnt;
                %fprintf(1,'StripNumber = %i, size(StripOffsets) = %i %i\n', TIF.StripNumber, size(TIF.StripOffsets));
            case 277         % sample_per pixel
                TIF.SamplesPerPixel  = entry.val;
%                 error( 'Color images not supported!' );
                %fprintf(1,'Color image: sample_per_pixel=%i\n',  TIF.SamplesPerPixel);
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
                %planar configuration is not fully supported here
%                 if ( TIF.PlanarConfiguration == 1 )
%                     error(sprintf('PlanarConfiguration = %i not supported', TIF.PlanarConfiguration));
%                 end
            case 296         % resolution unit
                IMG.resolution_unit= entry.val;
            case 305         % software
                IMG.software       = entry.val;
            case 306         % datetime
                IMG.datetime       = entry.val;
            case 315
                IMG.artist         = entry.val;
            case 317        %predictor for compression
                if (entry.val ~= 1), error('unsuported predictor value'); end
            case 320         % color map
                IMG.cmap           = entry.val;
                IMG.colors         = entry.cnt/3;
            case 339
                TIF.SampleFormat   = entry.val;
                if ( TIF.SampleFormat > 2 )
                    error(sprintf('unsuported sample format = %i', TIF.SampleFormat));
                end
            case 33628       %metamorph specific data
                IMG.MM_private1    = entry.val;
            case 33629       %this tag identify the image as a Metamorph stack!
                TIF.MM_stack       = entry.val;
                TIF.MM_stackCnt      = entry.cnt;
            case 33630       %metamorph stack data: wavelength
                TIF.MM_wavelength  = entry.val;
            case 33631       %metamorph stack data: gain/background?
                TIF.MM_private2    = entry.val;
            case 34412       % Zeiss LSM data (I have no idea what that represents...)
                % Jan: Image format of Zeiss Laser Scanning Microscope 510
                % Jan: apparently contains z-planes in low and full resolution in alternating order: full, low, full, low ...
                % Jan: do not assign this field to IMG since stack does not have this field, so "stack( img_read ) = IMG;" will not work
                % Jan: I don't know the meaning of the values in this LSM field anyway (it's a vector of uint8 numbers with length maxres)
                %IMG.LSM            = entry.val;
            otherwise
                fprintf(1,'ignored TIFF entry with tag %i (cnt %i)\n', TIF.entry_tag, entry.cnt);
        end
        % move to next IFD entry in the file
        fseek(TIF.file, file_pos+12,-1);
    end

    %total number of bytes per image:
    PlaneBytesCnt = IMG.width * IMG.height * TIF.BytesPerSample;

    if consolidateStrips
        %Try to consolidate the strips into a single one to speed-up reading:
        BytesCnt = TIF.StripByteCounts(1);

        if BytesCnt < PlaneBytesCnt

            ConsolidateCnt = 1;
            %Count how many Strip are needed to produce a plane
            while TIF.StripOffsets(1) + BytesCnt == TIF.StripOffsets(ConsolidateCnt+1)
                ConsolidateCnt = ConsolidateCnt + 1;
                BytesCnt = BytesCnt + TIF.StripByteCounts(ConsolidateCnt);
                if ( BytesCnt >= PlaneBytesCnt ), break; end
            end

            %Consolidate the Strips
            if ( BytesCnt <= PlaneBytesCnt(1) ) && ( ConsolidateCnt > 1 )
                %fprintf(1,'Consolidating %i stripes out of %i', ConsolidateCnt, TIF.StripNumber);
                TIF.StripByteCounts = [BytesCnt; TIF.StripByteCounts(ConsolidateCnt+1:TIF.StripNumber ) ];
                TIF.StripOffsets = TIF.StripOffsets( [1 , ConsolidateCnt+1:TIF.StripNumber] );
                TIF.StripNumber  = 1 + TIF.StripNumber - ConsolidateCnt;
            end
        end
    end

    %read the next IFD address:
    ifd_pos = fread(TIF.file, 1, 'uint32');

    if isfield( TIF, 'MM_stack' )

        if ( img_last > TIF.MM_stackCnt )
            img_last = TIF.MM_stackCnt;
        end
        
        %string description of the type of integer needed: int8 or int16...
        TIF.typecode = sprintf('int%i', TIF.BitsPerSample(1) );
        if ( TIF.SampleFormat == 1 )
            TIF.typecode = [ 'u', TIF.typecode ];
        end
        
        TIF.b = eval( [TIF.typecode '(0)'] );

        
        stack = zeros(Rwidth,Rheight, img_last-img_first,'uint16');
        %this loop is to read metamorph stacks:
        for ii = img_first:img_last

            TIF.StripCnt = 1;  %changed by read_plane!!

            %read the image
            fileOffset = PlaneBytesCnt * ( ii - 1 );

            img_read = img_read + 1;
            stack(:,:,img_read) = read_plane(fileOffset, IMG.width, IMG.height, 1);

            % UNCOMMENT THIS LINE TO GET MM STACK METADATA
            %[ IMG.info, IMG.MM_stack, IMG.MM_wavelength, IMG.MM_private2 ] = extractMetamorphData(ii);

        end

        break;  % stop reading directories

    else  %this part to read a normal TIFF stack:

        if ( img_skip + 1 >= img_first )

            TIF.StripCnt = 1;  %changed by read_plane!!
            %read the image

            % Jan: changes to read only the full resolution frames from lsm files
            if IMG.width == TIF.maxres
                
                img_read = img_read + 1;

                try
                    stack(:,:, img_read ) = read_plane(0, IMG.width, IMG.height, 1);
                catch
                    %stack
                    %IMG
                    error('The file contains dissimilar images: try to read them one by one');
                end
            end

        else
            img_skip = img_skip + 1;
        end

        if ( img_skip + img_read >= img_last )
            break;
        end
    end  %IF MM STACK
    
    
end %FOR EACH DIRECTORY ENTRY


%% distribute the MetaMorph info
if isfield(TIF, 'MM_stack') && isfield(TIF, 'info') && ~isempty(TIF.info)
    MM = parseMetamorphInfo(TIF.info, TIF.MM_stackCnt);
    exposure = [MM.Exposure];
end

% Build time axis
if exist('exposure','var'),
    nFrames = img_read;
    time = cumsum( [0 exposure(1:end-1)] );
else
    nFrames = img_read;
    time = 0:(nFrames-1);
end


%% clean-up
fclose(TIF.file);

drawnow;
%return empty array if nothing was read
if ~ exist( 'stack', 'var')
    stack = [];
end
clear IMG;
clear TIF;
return;


%============================================================================

function plane = read_plane(offset, width, height, planeCnt )
% WARNING: I have made the assumption that all the strips in each plane
% have been completely consolidated.  IE, each plane is a single strip.

global TIF;

%return an empty array if the sample format has zero bits
if ( TIF.BitsPerSample(planeCnt) == 0 )
    plane=[];
    return;
end

% Preallocate a matrix to hold the sample data:
%plane = repmat( TIF.b, width, height );

line = 1;

StripLength = TIF.StripByteCounts(TIF.StripCnt) ./ TIF.BytesPerSample(planeCnt);
assert( width*height==StripLength, 'Error: strips must be contiguous!' );

while ( TIF.StripCnt <= TIF.StripNumber )

    fseek(TIF.file, TIF.StripOffsets(TIF.StripCnt) + offset, 'bof');
    plane = fread( TIF.file, StripLength, TIF.typecode );  % );

    if ( length(plane) ~= StripLength )
        error('End of file reached unexpectedly.');
    end

    plane = reshape(plane, width, StripLength/width);
    
    TIF.StripCnt = TIF.StripCnt + 1;

    % copy the strip onto the data
    %plane(:, line:(line+size(strip,2)-1)) = strip;

    line = line + size(plane,2);
    if ( line > height )
        break;  %entire plane has been read, we're done
    end

end

% Extract valid part of data if needed
if ~all(size(plane) == [width height]),
  %  plane = plane(1:width, 1:height); NOC: 9-23-05
    error('Cropping data: more bytes read than needed...');
end

% transpose the image
plane = plane';

return;



%===================sub-functions that reads an IFD entry:===================


function [nbbytes, typechar] = matlab_type(tiff_typecode)
switch (tiff_typecode)
    case 1
        nbbytes=1;
        typechar='uint8';
    case 2
        nbbytes=1;
        typechar='uchar';
    case 3
        nbbytes=2;
        typechar='uint16';
    case 4
        nbbytes=4;
        typechar='uint32';
    case 5
        nbbytes=8;
        typechar='uint32';
    otherwise
        error('tiff type not supported')
end
return;

%===================sub-functions that reads an IFD entry:===================

function  entry = readIFDentry()

global TIF;
entry.typecode = fread(TIF.file, 1, 'uint16');
entry.cnt      = fread(TIF.file, 1, 'uint32');
%disp(['typecode =', num2str(entry.typecode),', cnt = ',num2str(entry.cnt)]);

[ entry.nbbytes, entry.typechar ] = matlab_type(entry.typecode);

if entry.nbbytes * entry.cnt > 4
    %next field contains an offset:
    offset = fread(TIF.file, 1, 'uint32');
    %disp(strcat('offset = ', num2str(offset)));
    fseek(TIF.file, offset, -1);
end

if TIF.entry_tag == 33629   %special metamorph 'rationals'
    entry.val = fread(TIF.file, 6*entry.cnt, entry.typechar);
else
    if entry.typecode == 5
        entry.val = fread(TIF.file, 2*entry.cnt, entry.typechar);
    else
        entry.val = fread(TIF.file, entry.cnt, entry.typechar);
    end
end
if ( entry.typecode == 2 ), entry.val = char(entry.val'); end

return;


%==============distribute the metamorph infos to each frame:
function [info, stack, wavelength, private2 ] = extractMetamorphData(imgCnt)

global TIF;

info = [];
stack = [];
wavelength = [];
private2 = [];

if TIF.MM_stackCnt == 1
    return;
end

left  = imgCnt - 1;

if isfield( TIF, 'info' )
    S = length(TIF.info) / TIF.MM_stackCnt;
    info = TIF.info( S*left+1:S*left+S );
end

if isfield( TIF, 'MM_stack' )
    S = length(TIF.MM_stack) / TIF.MM_stackCnt;
    stack = TIF.MM_stack( S*left+1:S*left+S );
end

if isfield( TIF, 'MM_wavelength' )
    S = length(TIF.MM_wavelength) / TIF.MM_stackCnt;
    wavelength = TIF.MM_wavelength( S*left+1:S*left+S );
end

if isfield( TIF, 'MM_private2' )
    S = length(TIF.MM_private2) / TIF.MM_stackCnt;
    private2 = TIF.MM_private2( S*left+1:S*left+S );
end


return;


%% %%  Parse the Metamorph camera info tag into respective fields
% EVBR 2/7/2005, FJN Dec. 2007
function mm = parseMetamorphInfo(info, cnt)

info   = regexprep(info, '\r\n|\o0', '\n');
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
                mm(k).Binning = sscanf(val, '%d x %d')';
            case 'Region'
                mm(k).Region = sscanf(val, '%d x %d, offset at (%d, %d)')';
            otherwise
                field  = regexprep(tok, ' ', '');
                if strcmp(val, 'Off')
                    eval(['mm(k).',field,'=0;']);
                elseif strcmp(val, 'On')
                    eval(['mm(k).',field,'=1;']);
                elseif isstrprop(val,'digit')
                    eval(['mm(k).',field,'=str2num(val)'';']);
                else
                    eval(['mm(k).',field,'=val;']);
                end
        end
    end
end

