function output = readTimecourse( filename )
% Reads data files from the fluorimeter, ignoring the header.
% Normally, there are many columns, one per experiment and each with its
% own time axis. The script takes the longest axis and assumes the
% integration time is the same for all. Injections are defined by lid
% openings, where the signal goes to zero. Each region of non-"zero"
% intensity is averaged for a datapoint in the output. Background ("zero")
% is separated from signal by the threshold defined below. This may need to
% be adjusted for each experiment.
%
% A list of files can also be chosen and these will be pooled together. The
% timecourses will be truncated to all be the same length. The time axes
% are assumed to be the same.
%

% Intensity threshold below which we assume the lid is open for adding the 
% titrated reagent. This separates the intensity data for each of the
% titration points.
% FIXME: may need to be changed for different experiments!
thresh = 10000;  % for uncorrected data!


% Get filenames from user
if nargin<1 || isempty(filename),
    [f,p] = uigetfile('*.txt','Select fluorescence timecourse file');
    if f==0, return; end
    filename = [p filesep f];
end



%% --- Read the header
fid = fopen( filename, 'r' );

% Read header for raw data.
d = textscan( fgetl(fid) ,'%f',1);
%nColumns = d{1}; %may be "<Trace>" instead... not useful.

d = textscan( fgetl(fid) ,'%f' );
nPoints = d{1};

d = textscan( fgetl(fid) ,'%s','Delimiter','\t');
header = d{1};

d = textscan( fgetl(fid) ,'%s');
fieldNames = d{1};
nColumns = numel(fieldNames);




%% --- Read the data
% We assume all x-axes are the same, but maybe different lengths.
d = textscan( fid, repmat('%f ',1,nColumns), 'Delimiter','\t' );
data = cell2mat(d);
% data = dlmread(filename);

fclose(fid);

% If there are multiple traces, select longest time axis column overall.
% Assume all files have the same step size (integration time).
[~,longestXaxis] = max( data(end,:) );
X = data(:,longestXaxis);

% Select data columns, skip other time axes.
Y = data(:,2:2:end);       

% Attempt to correct for photobleaching by dividing by the bleaching
% fraction over time. The bleaching curve was obtained by fitting the first
% 50 seconds of a number of curves to a line -- over this time there are no
% additions. The numbers vary, so I took an average. Here I assume the time
% axis is in seconds...
% TODO: this should be calculated from the data??
% TODO: for these figures, also draw big circles for the average titration
%   points detected and a plot of the result.
Yorig = Y;
Y = Y .* repmat(  (1+((0:size(Y,1)-1)*0.2273e-3))',  1, size(Y,2)  );  %comment out to disable bleaching correction

figure;
hUncorr = subplot(1,2,1);
plot(Yorig);
title('Uncorrected');
hCorr = subplot(1,2,2);
plot(Y); hold on;
title('Bleaching corrected');


% Find all points below a threshold of darkness and average each set of
% non-zero points to get the titration points.
output = zeros(1,size(Y,2));

for i=1:size(Y,2),
    % Seperate distinct events between the lid opennings.
    imagingRegions = Y(:,i)>thresh;
    [s,e] = rle( imagingRegions, 1 );
    
    % Skip the first and last datapoint of each span, which may include the
    % drop to and rise from zero when the lid is opened and closed.
    s = s+1; e = e-1;
    
    % Get the median of each span. Using the median helps avoid any spikes
    % or drops in intensity from opening/closing the lid or junk in the
    % cuvette.
    for j=1:numel(s),
        output(j,i) = median( Y(s(j):e(j),i) );
        
        % Plot on the raw fluorescence data to show the calculated mean value.
        scatter( (s(j)+e(j))/2.0, output(j,i), 50, 'ko','filled' );
    end
    
    % Normalize to initial fluorescence value
    output(:,i) = output(:,i)/output(1,i);
end
    

%% --- Save the output for importing into Origin.
outname = strrep(filename,'.txt','_analyzed.txt');
save( outname, 'output', '-ASCII' );






