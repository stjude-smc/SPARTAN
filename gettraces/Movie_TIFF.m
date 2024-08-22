classdef Movie_TIFF < Movie
% Movie_TIFF is a wrapper for a stack of TIFF images, as generated by
% imaging software like MetaMorph. See Movie for interface details.
% To load a movie, use the construtor with a filename. If a cell array of
% filenames is given, they are treated as a single movie that has been
% split with the first set of frames in the first file and the last set of
% frames in the last file. This is commonly done with very large files over
% 2 GB in size, which will not fit in the standard TIFF format.
%
%         movie = Movie_TIFF( filenames )
% 
% The readFrame(idx) and readFrames(indexes) methods can be used to read
% data from the movie file.

%   Copyright 2007-2024 All Rights Reserved.

% FIXME: does not correctly read XML-format metadata from MetaMorph.
% FIXME: gracefully handle DateTime field not found.
% 


%% ============ PROPERTIES ============= %
properties (SetAccess=protected, GetAccess=protected)
    swap = false;       %swap byte order if not the same as native (assumed to be little-endian)
end



%% ============ PUBLIC METHODS ============= %
methods
    
    % CONSTRUCTOR : path to .tif movie file
    function obj = Movie_TIFF( filename )

        % FIXME: need to call Movie_TIFF_MultiFile...
        if iscell(filename)
            filename = filename{1};
        end

        % Load TIFF metadata and save key tags.
        info = imfinfo( filename );
        obj.filename = filename;
        obj.header = info;
        obj.nFrames = numel(info);
        obj.nX = info(1).Width;
        obj.nY = info(1).Height;
        obj.precision = sprintf('uint%d',info(1).BitsPerSample);
        obj.swap = strcmpi(info(1).ByteOrder, 'big-endian');
        
        assert( info(1).SamplesPerPixel==1, 'Processing of multiple channels per pixel is not implemented' );
        
        
        % TIFF stacks from MicroManager/ImageJ may have only one IFD 
        % pointing to the entire image stack as one contiguous chunk.
        % These files contain two ImageDescription tags, but MATLAB can
        % only see one (starting with ImageJ=).
        if numel(info)==1 && numel(info.StripOffsets)==1 && isfield(info,'StripByteCounts')
            nFramesImplied = floor( (info.FileSize-info.StripOffsets)/info.StripByteCounts );
            if nFramesImplied>1
                obj.nFrames = nFramesImplied;
                obj.offsets = info.StripOffsets + (0:1:nFramesImplied)*info.StripByteCounts;
            end
        
        % If all frames are consolidated and uncompressed, can use fread directly.
        elseif numel(info(1).StripOffsets)==1 && strcmpi(info(1).Compression,'Uncompressed')
            obj.offsets = [info.StripOffsets];

        % Fall back on imread (slower but more flexible) 
        else
            obj.offsets = [];
        end


        % Determine actual measurement time axis:
        try
            obj.timeAxis = [];
            [ijExposure,ijText] = parseIJ(info);

            if isfield( info,'ExposureTime' )
                % EXIF tag 33434, used by FlashGordon
                ms = info(1).ExposureTime*1000;
                obj.timeAxis = (0:numel(info)-1)*ms;

            elseif ijExposure ~=0
                % ImageJ-specific metadata encoded in private tag 50839.
                obj.timeAxis = ijExposure*(0:obj.nFrames-1);
                obj.header.ImageJ = ijText;

            elseif isfield( info,'DateTime' )
                % Used by MetaMorph (others?)
                % NOTE: for most acquisition implementations, these
                % timestamps are for the time transfered from the camera
                % to the computer and may not actually be that valuable.
                dot = strfind( info(1).DateTime, '.' );
                if ~isempty(dot)
                    % This is hacky because ms time is variable width with
                    % MetaMorph and because each call to datenum is slow.
                    tcell = {info.DateTime};
                    ms  = cellfun( @(s) str2double(s(dot+1:end)), tcell );
                    pre = cellfun( @(s) s(1:dot-1), tcell, 'UniformOutput',false );
                    pre = vertcat( pre{:} );
                    times = ms + 24*60*60*1000*datenum(pre,'yyyymmdd HH:MM:SS')';

                    timediff = diff(times);
                    dt = median(timediff); %time resolution in ms
                    obj.timeAxis = dt*(0:obj.nFrames-1);
                end
            end
        catch
            obj.timeAxis = [];
        end
        
        % If no time information found, use frame numbers instead.
        % User will be asked in gettraces for a time resolution.
        if isempty( obj.timeAxis )
            disp('Warning: No exposure time found. Using frame numbers as time axis instead.');
            obj.timeAxis = 1:obj.nFrames;
        end
        
        % Get MetaMorph metadata. This can (in the STK format) be a special
        % tag with a whitespace delimited list of key/value pairs.
        if numel(info)==1 && isfield(info,'UnknownTags') && ...
                                      any( [info.UnknownTags.ID]==33628 )
            obj.header.MM = parseMetamorphInfo( info.ImageDescription, 1);
        end
        
        % Extract metadata from FlashGordon that specifies how to
        % subdivide image data into spectral channels.
        % OME-TIFF data will also be in this tag.
        if isfield(info,'ImageDescription') && startsWith(info(1).ImageDescription,'FlashGordon=')
            obj.metadata = parseFGImageDescription(info(1).ImageDescription);
        end
        
    end %constructor
    
    
    function data = readFrames( obj, idx )
        %assert( numel(idx)==1 && idx>=1 && idx<=obj.nFrames, 'Invalid index' ); 
        
        % Recursively get multiple frames if requested.
        if numel(idx)>1
            data = zeros( obj.nY,obj.nX,numel(idx), obj.precision );

            for i=1:numel(idx)
                data(:,:,i) = obj.readFrames( idx(i) );
            end
            return;
        end
        
        if ~isempty(obj.offsets)
            % Fast version for optimized TIFF stacks.
            fid = fopen( obj.filename, 'r' );
            fseek( fid, obj.offsets(idx), -1 );
            data = fread( fid, [obj.nX obj.nY], ['*' obj.precision] )';
            if obj.swap, data = swapbytes(data); end
            fclose(fid);
        else
            % Slower version tolerant to fragmentation and unusual parameters.
            data = imread( obj.filename, 'Info',obj.header, 'Index',idx );
        end
    end
    
end %public methods



end %class Movie_TIFF






%% ============ ACCESSORY FUNCTIONS ============= %

function [Interval_ms,text] = parseIJ(info)
% Parses TIFF tags 50838 and 50839 associated with MicroManager/ImageJ.
% INPUT: TIFF metadata struct from imfinfo.
% OUTPUT: standardized metadata struct (see next function).
% See: https://micro-manager.org/Micro-Manager_File_Formats
% Tag 50839 includes JSON text after a brief binary header.
% The details are unclear, so I look for { that starts the JSON.
% For now, we only grab the time resolution if available

Interval_ms = 0;

try
    for i=1:numel(info.UnknownTags)
        val = info.UnknownTags(i).Value;
    
        if info.UnknownTags(i).ID == 50839
            if ~strcmp( char(val(1:4)), 'IJIJ' ), return; end
            idxBegin = find(val == '{',1,'first');
            if isempty(idxBegin), return; end
            text = char(val(idxBegin:2:end));
            data = jsondecode( text );
            Interval_ms = data.Interval_ms;
        end
    end
catch
    Interval_ms = 0;
end

end %function parseIJ



function metadata = parseFGImageDescription(input)
% Image Description (tag 270) text is saved by Flash Gordon to assist in
% identifying spectral channels that are tiled together in each image.

    metadata = [];
    output = [];
    validFields = {'FlashGordon','hardware','binning','exposureTime','frameTime','StageX','StageY','fieldArrangement','power_mW'};
    fieldFormat = {'%s','%s','%s','%f','%f','%f','%f',[],'%f'};
    
    validChFields = {'name','wavelength','photonsPerCount'};
    chFieldFormat = {'%s','%d','%f'};
    
    validLaserFields = {'wavelength','dutyCycle','framesActive'};
    laserFieldFormat = {'%d','%f','%s'};
    
    % Read lines as name=value pairs, skipping unrecognized fields.
    input = splitlines(input);
    for i=1:numel(input)
        line = split( input{i}, '=' );
        if numel(line)~=2, continue; end
        [fieldName,val] = line{:};
        
        try
            % Build channels struct array
            if startsWith(fieldName,'channel') && fieldName(9)=='.'
                chID = str2double(fieldName(8));
                subfield = fieldName(10:end);
                formatID = strcmp(subfield,validChFields);
                output.channels(chID).(subfield) = sscanf( val, chFieldFormat{formatID} );
                
            % Illumination struct array
            elseif startsWith(fieldName,'laser') && fieldName(7)=='.'
                laserID = str2double(fieldName(6));
                subfield = fieldName(8:end);
                formatID = strcmp(subfield,validLaserFields);
                
                if strcmpi(subfield,'framesActive')
                    assert( all(ismember(val,'1234567890[],;: ')), 'Invalid framesActive' );
                    output.lasers(laserID).(subfield) = eval(val);
                else
                    output.lasers(laserID).(subfield) = sscanf( val, laserFieldFormat{formatID} );
                end
                

            % FieldArrangement is text to define a matlab-style matrix
            % that describes how to split image data into channels.
            elseif strcmpi(fieldName,'fieldArrangement')
                assert( all(ismember(val,'1234567890[],;: ')), 'Invalid field arrangement' )
                output.fieldArrangement = eval(val);
            else
                formatID = strcmp(fieldName,validFields);
                output.(fieldName) = sscanf( val, fieldFormat{formatID} );
            end
        catch
            disp(['Ignoring invalid FlashGordon ImageDescription line: ' input{i}]);
        end
    end
    
    % Verify channel information is valid so it can be used by ChannelExtractor.
    try
        nCh = numel(output.channels);
        fa = output.fieldArrangement(:);
        wavelength = [output.channels.wavelength];
        ppc = [output.channels.photonsPerCount];
        
        if ~all(fa>=0 & fa<=nCh & sum(fa>0)==nCh)            || ...
           any(isnan(wavelength)) || numel(wavelength)~=nCh  || ...
           any(isnan(ppc)) || numel(ppc)~=nCh                || ...
           any(cellfun(@isempty,{output.channels.name}))
           error('Invalid channel metadata');
        end
    catch
        warning('FlashGordon metadata not valid; ignoring');
        return;
    end
    
    metadata = output;  %only return a value if we pass all tests.
    
end %function




% Parse the Metamorph camera info tag into respective fields
% This is only used with old-format (STK) movies from MetaMorph.
% EVBR 2/7/2005, FJN Dec. 2007. See tiffread.m
function mm = parseMetamorphInfo(info, cnt)

% info   = regexprep(info, '\r\n|\o0', '\n');
info( info==char(0) | info==char(13) ) = newline;

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
        [v, ~, ~, pos] = sscanf(val, '%i');
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
