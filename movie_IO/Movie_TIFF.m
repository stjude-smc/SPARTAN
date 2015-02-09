classdef Movie_TIFF
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
%

% FIXME: does not correctly read XML-format metadata from MetaMorph.
% FIXME: gracefully handle DateTime field not found.
% 


%% ============ PROPERTIES ============= %
properties (SetAccess=protected, GetAccess=public)
    filenames = {}; % full path and filename to loaded file
    filetag = '';   % original base filename of movie file series
    
    nX=0;       % size (in pixels) of x dimension (columns)
    nY=0;       % size (in pixels) of y dimension (rows)
    nFrames=0;  % number of images in the stack
    
    precision='';  %pixel data format and bitrate (usually uint16)
    
    timeAxis=[];     % uniform time axis in ms, relative to start.
    timestamps = []; % actual movie timestamps in ms, relative to start.
    
    header = struct([]);  % TIFF metadata for first frame.

end %end public properties


properties (SetAccess=protected, GetAccess=protected),
    movieHeaders = [];  %metadata for all frames across all files (struct).
    nFramesPerMovie = [];
    
    offsets = {};  %for each file, list of byte offsets to binary frame data.
end



%% ============ PUBLIC METHODS ============= %
methods
    
    % Basic constructor with the filename of a TIFF movie file to load or a
    % cell array of filenames if the movie was split into multiple files.
    % This may happenen with very large movies (>2GB).
    function obj = Movie_TIFF( filenames )
        
        if ~iscell(filenames),  filenames = {filenames};  end
        obj.filenames = filenames;
        nFiles = numel(filenames);
        
        % Extract TIFF tags (metadata) from all files.
        % This will be used to help determine the order.
        headers = cell(1,nFiles);
        times = cell(nFiles,1);  %time axis for each movie.
        firstTime = zeros(nFiles,1); %time at first frame in each file.
        
        for i=1:nFiles,
            % Get TIFF tag metadata for this file.
            info = imfinfo( filenames{i} );
            headers{i} = info;
            obj.nFramesPerMovie(i) = numel(info);
            
            if numel(info(1).StripOffsets)==1
                obj.offsets{i} = [info.StripOffsets];
            end
            
            % FIXME: verify all required field names are present.
            assert( info(1).SamplesPerPixel==1, 'Processing of multiple channels per pixel is not implemented' );
            
            % Process time axis information.
            if isfield( info,'DateTime' )
                dot = strfind( info(1).DateTime, '.' );
                
                % NOTE: for most acquisition implementations, these
                % timestamps are for the time transfered from the camera 
                % to the computer and may not actually be that valuable.
                % This is hacky because ms time is variable width with
                % MetaMorph and because each call to datenum is slow.
                tcell = {info.DateTime};
                ms  = cellfun( @(s) str2double(s(dot+1:end)), tcell );
                pre = cellfun( @(s) s(1:dot-1), tcell, 'UniformOutput',false );
                pre = vertcat( pre{:} );    
                times{i} = ms + 24*60*60*1000*datenum(pre,'yyyymmdd HH:MM:SS')';

                firstTime(i) = times{i}(1);
                
            elseif isfield( info,'ExposureTime' )
                % No explicit time axis is available, but we can
                % reconstruct an approximate one using the exposure time.
                % This field is specific to the LabVIEW software.
                ms = info(1).ExposureTime*1000;                
                firstTime(i) = 0 + sum( obj.nFramesPerMovie(1:i-1) )*ms;
                times{i} = firstTime(i) + (0:numel(info)-1)*ms;
                
            else
                warning('Movie_TIFF:NoDateTime','No DateTime field. Using frame number as time axis.');
                firstTime(i) = 1 + sum( obj.nFramesPerMovie(1:i-1) );
                times{i} = firstTime(i)+(1:numel(info))-1;
            end
            
            % Parse file dimensions, etc
            if i==1,
                [p,f] = fileparts(info(1).Filename);
                f = regexprep(f,'-file[0-9]*$','');
                obj.filetag = fullfile(p,f);
                
                obj.nX = info(1).Width;
                obj.nY = info(1).Height;
                
                obj.precision = sprintf('uint%d',info(1).BitsPerSample);
            else
                % Verify dimensions are the same
                assert( info(1).Width==obj.nX && info(1).Height==obj.nY, ...
                           'Movies in series have different sizes!' );
                 
                % Verify the file base name is the same across files.
                % Remove the extension and the '-file003' prefix.
                % We expect (but don't assume) the extension is '.tif'.
                [p,f] = fileparts(info(1).Filename);
                f = regexprep(f,'-file[0-9]*$','');
                tag = fullfile(p,f);
                
                if ~strcmp(tag,obj.filetag),
                   warning('Movie_TIFF:filename_mismatch', ...
                        'Files do not appear to be part of the same acquisition');
                    disp(info(1).Filename);
                end
            end
        end
        
        % Reorder the files so the frames are continuous in time.
        % This is not generally necessary because sorting the filenames
        % gives the correct order.
        assert( numel(unique(firstTime)) == numel(firstTime), ...
                'Not enough information to put movie files in order!' );
        
        [~,fileOrder] = sort(firstTime);
        obj.filenames       = obj.filenames(fileOrder);
        obj.nFramesPerMovie = obj.nFramesPerMovie(fileOrder);
        obj.movieHeaders    = headers(fileOrder); %convert to matrix
        obj.timestamps      = [ times{fileOrder}   ]; %convert to matrix
        obj.header = obj.movieHeaders{1}(1);
        
        % Construct a time axis with a uniform time resolution (.timeAxis).
        % Some programs will be confused by the actual timestamps with
        % non-uniform time resolution.
        obj.nFrames  = sum( obj.nFramesPerMovie );
        timediff = diff(obj.timestamps);
        
        if obj.timestamps(1)==1 && all(timediff==1),
            % Time in frames
            obj.timeAxis = obj.timestamps;
        else
            dt = round( 10*median(timediff) )/10; %time resolution in ms
            obj.timeAxis = 0:dt:(dt*obj.nFrames-1);
            
            % Verify the timestamps are continuous. If files from distinct
            % movies are spliced together, they will have big jumps. This is a
            % warning because sometimes movies are shuttered, giving big jumps.
            if any(timediff<0.1) || any(timediff>3*dt),
                warning('Movie_TIFF:discontinuousTimestamps', ...
                        'Timestamps are not continuous. This can happen if files from multiple movies are mixed!');
                disp(  [ min(diff(obj.timestamps)) max(diff(obj.timestamps)) ]  );
            end
            %obj.timestamps = obj.timestamps-obj.timestamps(1); %relative to start.
        end
        
                
        % Get MetaMorph metadata. This can (in the STK format) be a special
        % tag with a whitespace delimited list of key/value pairs.
        if numel(info)==1 && isfield(info,'UnknownTags') && ...
                                      any( [info.UnknownTags.ID]==33628 ),
                                  
            obj.header.MM = parseMetamorphInfo( info.ImageDescription, 1);
        end
        
        % In new versions (TIFF), the metadata is instead in XML.
        if isfield(obj,'ImageDescription'),
            % TODO
        end
        
        % If using fread for optimized files, no need to keep headers. This can
        % speed up the parfor loop in gettraces by minimizing data transfers.
        if ~isempty(obj.offsets) && all( ~cellfun(@isempty,obj.offsets) ),
            obj.movieHeaders = {};
        end
        
    end %constructor
    
    
    
    
    % Data access methods. Data are only loaded when needed by these functions.
    function data = readFrames( obj, idx )
        
        % Parse input arguments; insure correct orientation of vector.
        idx = reshape(idx, [1 numel(idx)]);
        assert( min(idx)>=1 && max(idx)<=obj.nFrames, 'Invalid indexes' );
        
        % Preallocate space for the output data. (Removing star in precision)
        data = zeros( obj.nY,obj.nX,numel(idx), obj.precision );
        
        for i=1:numel(idx),
            data(:,:,i) = obj.readFrame( idx(i) );
        end
    end
    
    
    function data = readFrame( obj, idx )
        assert( numel(idx)==1 && idx>=1 && idx<=obj.nFrames, 'Invalid index' ); 
        
        % Determine which file this frame number belongs to.
        movieFirstFrame = 1+cumsum([0 obj.nFramesPerMovie]);
        idxFile = find( idx>=movieFirstFrame, 1, 'last' );
        idx = idx - movieFirstFrame(idxFile)+1;
        
        %offsets = obj.movieHeaders{idxFile}(idx).StripOffsets;
        
        if isempty(obj.offsets) || isempty(obj.offsets{idxFile}),
            % Slower version that can read frames broken up into strips.
            data = imread( obj.filenames{idxFile}, ...
                     'Info',obj.movieHeaders{idxFile}, 'Index',idx );
        else
            % Fast version for optimized TIFF stacks.
            fid = fopen( obj.filenames{idxFile}, 'r' );
            fseek( fid, obj.offsets{idxFile}(idx), -1 );
            data = fread( fid, [obj.nX obj.nY], ['*' obj.precision] )';
            fclose(fid);
        end
    end
    
end %public methods



end %class Movie_TIFF






%% ============ ACCESSORY FUNCTIONS ============= %

% Parse the Metamorph camera info tag into respective fields
% This is only used with old-format (STK) movies from MetaMorph.
% EVBR 2/7/2005, FJN Dec. 2007. See tiffread.m
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
