function [stkData,peaks,image_t] = gettraces(varargin)
% GETTRACES  Extract smFluorescence traces from movies
%
%    GETTRACES()
%    Launches the gettraces graphical user interface.
%
%    [STK,PEAKS,IMG] = GETTRACES( FILENAME, PARAMS )
%    Loads a movie from FILENAME, finds fluorescence peaks, and returns
%    their locations (PEAKS) as a Nx2 matrix (x,y).  IMG is the image used
%    for selecting intensity peaks.  PARAMS (optional) specifies the
%    intensity threshold (don_thresh) and minimal peak seperation
%    (overlap_thresh).
%
%    [STK,PEAKS,IMG] = gettraces( FILENAME, PARAMS, OUTFILE )
%    As above, but also extracts traces from the PEAKS locations and saves
%    the resulting traces to OUTFILE.
%
%    STK can also be passed instead of FILENAME if the movie has already
%    been loaded by gettraces. This saves on load time.
%
%    GETTRACES( DIRECTORY, PARAMS )
%    For each filename in DIRECTORY, the movie is loaded, processed, and
%    a .traces file is saved automatically. No data is returned.
%    This is "batch mode".  Additional (optional) fields in PARAMS include
%    option to check all child folders (recursive) and option to skip
%    movies that have already been processed (skipExisting).
%

% TODO: also return unfiltered peak list (overlap=0) so overlap statistics
% can be displayed in scripts that call this function.


%% Parse input arguments

% If calling directly from command line,
% launch the GUI.
if nargin==0,
    gettraces_gui;
    return;
end



    
%------ If a structure is specified, load as STK data
if nargin>=1 && isstruct(varargin{1}),
    stkData = varargin{1};

%------ Otherwise, load an individual file
elseif nargin>=1 && ischar(varargin{1}) && ~isdir(varargin{1}),
    stkData = OpenStk2( varargin{1} );

%------ If a directory name is given, run in batch mode
elseif nargin>=1 && ischar(varargin{1}) && isdir(varargin{1})
    if nargin>=2,
        params = varargin{2};
        batchmode( varargin{1}, params );
    else
        batchmode( varargin{1} );
    end
    
    return;

else
    error('gettraces: Invalid param 1');
end




%------ Find peaks of intensity corrosponding to isolated single molecules
params.don_thresh = 0;
params.overlap_thresh = 2.1;

if nargin>=2,
    params = varargin{2};
end
    
% Generate image that will be used to select peaks
image_t = stkData.stk_top - stkData.background;
[nrow ncol] = size(image_t);

if ~isfield(params,'don_thresh') || params.don_thresh==0
    if ~isfield(params,'thresh_std')
        thresh_std = 8;
    else
        thresh_std = params.thresh_std;
    end

    % Automatically threshold setting based on variance of
    % background intensity at end of movie.
    endBG = sort( stkData.endBackground(:) );
    endBG_lowerHalf = endBG( 1:floor(numel(endBG)*0.75) );
    don_thresh = thresh_std*std( endBG_lowerHalf );
else
   don_thresh = params.don_thresh-mean2(stkData.background);
end

overlap_thresh = 2.1;
if isfield(params,'overlap_thresh')
    overlap_thresh = params.overlap_thresh;
end

% Find peak locations from total intensity
[peaksX,peaksY] = getPeaks( image_t, don_thresh, overlap_thresh );
peaks = [peaksX peaksY];



%------ Extract traces using peak locations
if nargin>=3 && ischar(varargin{3}),
    outputFilename = varargin{3};
    integrateAndSave( stkData.stk, stkData.stk_top, peaks', ...
        outputFilename, stkData.time );
end
    
return;


%% --------------------- OPEN STK MOVIE (CALLBACK) --------------------- %



function [stkData] = OpenStk2(filename)

% If the movie is compressed, deflate it to a temporary location first
if strfind(filename,'.stk.bz2'),
    
    % decompress
    z_fname = filename;
    status = system( ['bunzip2 --keep "' z_fname '"'] );
    if status ~= 0, 
        error('Error uncompressing STK: %s',z_fname);
    end
    
    % Load stk movie -- stk(X,Y,frame number)
    filename = strrep(z_fname,'.stk.bz2','.stk');
    [stk,time] = tiffread(filename);
    [stkX,stkY,Nframes] = size(stk);
    
    % Delete the compressed movie -- we no longer need it!
    delete( z_fname );

elseif strfind(filename,'.movie'),
    fid = fopen( filename, 'r' );
    stkSize = fread( fid, 3, 'int16' );
    stkX = stkSize(1); stkY = stkSize(2); Nframes = stkSize(3);
%     time = fread( fid, Nframes, 'float' );
    stk = fread( fid, stkX*stkY*Nframes, 'int16' );
    stk = reshape(stk,stkSize');
    fclose(fid);
    time = 1:Nframes;
    
else
    % Load stk movie -- stk(X,Y,frame number)
    [stk,time] = tiffread(filename);
    [stkX,stkY,Nframes] = size(stk);
end


% Create an average image of the first 10 frames (stk_top)
averagesize = min([10 Nframes]);
movie = stk(:,:,1:averagesize);
stk_top = mean(movie,3);



% Create an estimated background image by:
% 1. Divide the image into den*den squares
% 2. For each square, find the fifth lowest number
% 3. Rescaling these values back to the original image size
if stkX == 128
    den=4;  %not optimized
elseif stkX == 170
    den=5;
%     set(handles.overlap,'String', 2.1);
elseif stkX == 256
    den=8;
%     set(handles.overlap,'String', 2.5);
elseif stkX == 512
    den=16;
%     set(handles.overlap,'String', 4.5);
end

background = stk_top;  %**
temp = zeros( stkY/den, stkX/den );

for i=1:stkY/den
    for j=1:stkX/den
        sort_temp = background(den*(i-1)+1:den*i,den*(j-1)+1:den*j);
        sort_temp = reshape(sort_temp,1,den*den);  % make into a vector
        sort_temp = sort(sort_temp);

        temp(i,j) = sort_temp(den);  % get the 1/den % smallest value
    end
end

% Rescale the image back to the actual size
background=imresize(temp,[stkX stkY],'bilinear');
% handles.background = mean( stk(:,:,end-4:end), 3 );

% Also create a background image from the last few frames;
% useful for defining a threshold.
endBackground = double(  stk(:,1:stkX/2,       end-11:end-1)   ...
                      +  stk(:,(1+stkX/2):end, end-11:end-1)   );

stkData.stk = stk;
stkData.stk_top = stk_top;
stkData.background = background;
stkData.endBackground = endBackground;
stkData.time = time;

[stkData.stkX,stkData.stkX,stkData.stkX] = size(stk);
stkData.nFrames = size(stk,3);

% END FUNCTION OpenStk




% --------------- OPEN ALL STK'S IN DIRECTORY (BATCH) --------------- %

function batchmode(direct,params)

if nargin<2,
    params = struct([]);
end

if ~isfield(params,'skipExisting'),
    params(1).skipExisting = 0;
end
if ~isfield(params,'recursive')
    params.recursive = 0;
end
if ~isfield(params,'don_thresh'),
    params.don_thresh = 0;
end
if ~isfield(params,'overlap_thresh'),
    params.overlap_thresh = 2.1;
end

% Create header for log file
log_fid = fopen( [direct filesep 'gettraces.log'], 'w' );
fprintf(log_fid, 'Donor Thresh = %.1f\n',params.don_thresh);
fprintf(log_fid, 'Overlap = %.1f\n',params.overlap_thresh);
% fprintf(log_fid, 'Acceptor Thresh = %.1f\n',acc_thresh);

fprintf(log_fid,'\n%s\n\n%s\n%s\n\n%s\n',date,'DIRECTORY',direct,'FILES');

% Get list of files in current directory (option: and all subdirectories)
if params.recursive
    stk_files  = rdir([direct filesep '**' filesep '*.stk*']);
else
    stk_files  = rdir([direct filesep '*.stk*']);
end


% For each file in the user-selected directory
i = 0;
nFiles = length(stk_files);
for file = stk_files'
    
    i = i+1;
    
    % Skip if previously processed (.traces file exists)
    stk_fname = strrep(file.name,'.bz2','');
    [p,name]=fileparts(stk_fname);
    traceFname = [p filesep name '.traces'];
    
    if params.skipExisting && exist(traceFname,'file'),
        %disp( ['Skipping (already processed): ' stk_fname] );
        fprintf(log_fid, 'Skip %s\n', file.name);
        continue;
    end
    
    if ~exist('h','var')
        h = waitbar(0,'Extracting traces from movies...');
    end
    
    % Load STK file
    stkData = OpenStk2( file.name );
    
    % Pick molecules using default parameter values
    image_t = stkData.stk_top - stkData.background;
    [nrow ncol] = size(image_t);

    if ~isfield(params,'don_thresh') || params.don_thresh==0
        if ~isfield(params,'thresh_std')
            thresh_std = 8;
        else
            thresh_std = params.thresh_std;
        end
        
        % Automatically threshold setting based on variance of
        % background intensity at end of movie.
        endBG = sort( stkData.endBackground(:) );
        endBG_lowerHalf = endBG( 1:floor(numel(endBG)*0.75) );
        don_thresh = thresh_std*std( endBG_lowerHalf );
    else
        don_thresh = params.don_thresh-mean2(stkData.background);
    end

    % Find peak locations from total intensity
    [peaksX,peaksY] = getPeaks( image_t, don_thresh, params.overlap_thresh );
    peaks = [peaksX peaksY];
    
    % Save the traces to file
    integrateAndSave( stkData.stk, stkData.stk_top, peaks', ...
        traceFname, stkData.time );
    
    % Save entry in log file
    fprintf(log_fid, '%.0f %s\n', size(peaks,1)/2, file.name);
    
    waitbar(i/nFiles); drawnow;
end

if exist('h','var')
    close(h);
end

fclose(log_fid);







% --------------- PICK MOLECULES CALLBACKS --------------- %

%------------- Pick Cy3 spots ----------------- 
function [picksX,picksY] = getPeaks( image_t, threshold, overlap_thresh )
% Localizes the peaks of molecules in the Cy3 channel and infers the
% positions of the Cy5 peaks by applying a transofmration to the Cy3 peak
% positions.
%
% picksX - X-coords of all molcules, as Cy3,Cy5,Cy3,Cy5 in order
% picksY - Y-coords ...

[nrow ncol] = size(image_t);

%---- 1. Get coordinates of donor-side peaks
donor_t = image_t(:,1:ncol/2);
acceptor_t = image_t(:,(ncol/2)+1:end);
total_t = donor_t+acceptor_t;

[don_x,don_y] = pickPeaks( total_t, threshold, overlap_thresh );
nPicked = numel(don_x);

acc_x = don_x + (ncol/2);
acc_y = don_y;


%---- 6. Output results

picksX = zeros(nPicked*2,1);
picksY = zeros(nPicked*2,1);

for i=1:nPicked,
    picksX(2*i-1) = don_x(i);
    picksY(2*i-1) = don_y(i);
    picksX(2*i)   = acc_x(i);
    picksY(2*i)   = acc_y(i);
end    

return;  %skip everything else!


%---- 2. Estimate coordinates of acceptor-side peaks

% Initial guess: simple translation to right half (no rotation, scale)
acc_x = don_x + (ncol/2);
acc_y = don_y;

for j=1:nPicked,

    % Refine acceptor peak positions by finding local maxima
    % within the 3x3 grid around the initial guess.
    temp = image_t( acc_y(j)-1:acc_y(j)+1, acc_x(j)-1:acc_x(j)+1 );
    [maxy, maxx] = find(temp==max(max(temp)),1);
    acc_x(j) = acc_x(j) +maxx-2;
    acc_y(j) = acc_y(j) +maxy-2;
    
end

% OR do the following:
% help says do not use if scale changes in transform
%acceptor_points = cpcorr( acceptor_points, donor_points, acceptor_t, donor_t );


%---- 3. Create mapping function: donor to acceptor side
% ...to account for misalignment of DualView.  The transform is then 
% applied to get predicted Cy5 peak positions.
% 
%   Tranform matrix format:
%   scale X | shear Y | 0
%   shear X | scale Y | 0
%   trans X | trans Y | 1
% 
% NOTE that the image is a grid, so transforms of non-integer amounts
% have stochastic effects going from continuous to discrete cooredinates.
% 
% AFFINE works best
% 

if nPicked<5, 
    % If less than 5 molecules, don't attempt to derive a transform,
    % just use translation
    disp( 'Warning: Too few molecules to create transformation' );
    f = [1 0 0 ; 0 1 0 ; ncol/2 0 1];
else
    donor_points    = [don_x ; don_y]';
    acceptor_points = [acc_x ; acc_y]';

    tform = cp2tform(acceptor_points,donor_points,'linear conformal');  %image proc TK
    f = tform.tdata.Tinv;
    
%     u = [0 1]; v=[0 0];
%     [tx,ty] = tformfwd(tform, u,v);
%     dx = tx(2)-tx(1)
%     dy = ty(2)-ty(1)
%     angle = (180/pi) * atan2(dy, dx) 
%     scale = 1 / sqrt(dx^2 + dy^2)
end



%---- 4. Use the transform to get acceptor-side coordinates
for i=1:nPicked,
    
    % Apply transform to get acceptor-side coordinates
    points=[don_x(i) don_y(i) 1]*f;
    points=round(points);

    acc_x(i) = points(1);
    acc_y(i) = points(2);
end    


%---- 5. Repick based on total peak intensity
% Picking based on donor intensity alone may preferentially miss
% traces with high average FRET (low donor intensity)

% reg_acceptor = imtransform( acceptor_t, tform, 'Fill',255 );
% reg_acceptor = reg_acceptor(1:170,1:85);





%---- 6. Output results

picksX = zeros(nPicked*2,1);
picksY = zeros(nPicked*2,1);

for i=1:nPicked,
    picksX(2*i-1) = don_x(i);
    picksY(2*i-1) = don_y(i);
    picksX(2*i)   = acc_x(i);
    picksY(2*i)   = acc_y(i);
end    



% end function getTraces



function wavg = WeightedAvg( vals, weights )
% Produces a weighted average of <vals>.
% Used by getTraces to find fluor centroid position

weights = weights / sum(weights);  % normalize

wavg = sum( vals.*weights );

% end function WEIGHTEDAVG



function [don_x,don_y] = pickPeaks( image_t, threshold, overlap_thresh )
% Localizes the peaks of fluorescence, removing any too close together.
% NOTE: only the donor side is passed to this function

[nrow ncol] = size(image_t);


% 1. For each pixel (excluding edges), pick those above threshold that have
% greater intensity than their neighbors (3x3,local maxima) -- donor only.
% tempxy is peak position, centroidxy is estimated true molecule position.
nMols=0;
for i=1+3:nrow-3,
    for j=1+3:ncol-3,
        block = image_t(i-1:i+1,j-1:j+1);
        cross = block( [2,4,5,6,8] );
        
        if image_t(i,j)>threshold && image_t(i,j)==max(cross)
            
            % Calc centroid position using intensity-weighted average of
            % position
            off_x = WeightedAvg( 1:3, sum(block,1)  ) -1.5;
            off_y = WeightedAvg( 1:3, sum(block,2)' ) -1.5;
%                   
%             line(j+off_x-0.5,i+off_y-0.5,'marker','o','color','y','EraseMode','background')
            
            % Save the position of this molecule
            nMols=nMols+1;
            tempx(nMols)=j;
            tempy(nMols)=i;
            
            centroidx(nMols) = j+ off_x;
            centroidy(nMols) = i+ off_y;
        end
    end
end

if nMols<1,
    don_x = [];
    don_y = [];
    return;
end


% 2. Remove pick pairs that are closer than the cutoff radius:
% Signals may have cross-contamination. -- donor only
% THIS STEP NEEDS IMPROVEMENT - if molecules are too close, they will not
% be picked as seperate peaks.
overlap=zeros(1,nMols);

for i=1:nMols,
    for j=i+1:nMols,
        % quickly ignore molecules way beyond the threshold
        if abs(centroidy(j)-centroidy(i)) > overlap_thresh
            break;
        end
        % otherwise, we have to check more precisely
        if (centroidx(i)-centroidx(j))^2 + (centroidy(i)-centroidy(j))^2 ...
           < overlap_thresh^2
            overlap(i)=1;
            overlap(j)=1;
        end
    end
end


indexes = find( overlap==0 ); %find peaks with overlap less than threshold
don_x = tempx(indexes);
don_y = tempy(indexes);


% END FUNCTION pickPeaks






% --------------------- SAVE PICKED TRACES TO FILE --------------------- %


% function saveTraces( handles, hObject )
function integrateAndSave( stk, stk_top, peaks, stk_fname, time )
% NOTE: can find which pixels to use by correlation (original code did so?)

[stkX,stkY,Nframes] = size(stk);

if stkX == 128
    squarewidth=3;
    NumPixels=5;  %unoptimized
elseif stkX == 170
    squarewidth=3;
    NumPixels=4;
elseif stkX == 256
    squarewidth=3;
    NumPixels=4;
elseif stkX == 512
    squarewidth=5;
    NumPixels=16;
end

% Get x,y coordinates of picked peaks
Npeaks = size(peaks,2)/2;
x = peaks(1,:);
y = peaks(2,:);

regions=zeros(NumPixels,2,2*Npeaks);  %pixel#, dimension(x,y), peak#

% Define regions over which to integrate each peak --
% Done separately for each channel!
for m=1:2*Npeaks
    
    hw = floor(squarewidth/2);
    
    % Get pixels around picked point (squarewidth^2 pixels)
    peak = stk_top( ...
            y(m)-hw:y(m)+hw, ...
            x(m)-hw:x(m)+hw  );
    center = sort( peak(:) );  %vector of sorted intensities
    
    % Get pixels whose intensity is greater than the median (max=NumPixels).
    % We just want the centroid to avoid adding noise ...
    % A is x-coord, B is Y-coord of the top <NumPixels> pixels
    [A,B]=find( peak>=center(squarewidth*squarewidth-NumPixels+1), NumPixels );
    
    % Define a region over which to integrate each peak
    regions(:,:,m) = [ A+y(m)-hw-1, B+x(m)-hw-1  ];
end


% Create a trace for each molecule across the entire movie
traces = zeros(2*Npeaks,Nframes);

s = size(stk);
idx = sub2ind( s(1:2), regions(:,1,:), regions(:,2,:) );
% bg = handles.background;

for k=1:Nframes,
    frame = double(stk(:,:,k));
    traces(:,k) = sum( frame(idx) );
end

clear stk;


% Subtract background from last point of all traces so data
% can fit into an int16 matrix
donor    = traces(1:2:end,:);
acceptor = traces(2:2:end,:);

donor    = donor    - mean( donor(:,end)    );
acceptor = acceptor - mean( acceptor(:,end) );


% Save data to file.
stk_fname = strrep(stk_fname,'.bz2','');
[p,name]=fileparts(stk_fname);
save_fname = [p filesep name '.traces'];

saveTraces( save_fname, 'traces', donor,acceptor, [], time );


% Save the locations of the picked peaks for later lookup.
% FORMAT:  mol_name, don x, don y, acc x, acc y
% filename=strrep(handles.stkfile,'.stk','.loc.txt');
% fid = fopen(filename,'w');
% 
% for j=1:nTrances/2;
%     don_x = x(2*j-1);
%     don_y = y(2*j-1);
%     acc_x = x(2*j);
%     acc_y = y(2*j);
%     
%     fprintf(fid, '%s_%d- %d %d %d %d\n', name, j, don_x,don_y,acc_x,acc_y);
% end
% 
% fclose(fid);


