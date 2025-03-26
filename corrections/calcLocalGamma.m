function result = calcLocalGamma(files)
% calcLocalGamma  spatial map of local variance in detection efficiency.
%
%   R = calcLocalGamma(FILES) calculates a spatial mapping function R(x,y)
%   of the acceptor channel detection sensitivity relative to the donor
%   channel, which can be used by applyLocalGamma to correct for
%   non-uniform detection efficiency (flat-field correction).
%   
%   FILES can be a char array path to a .traces file, or cell array of 
%   such paths, obtained from movies with a standard sample, such as DNA
%   oligos, with a fixed FRET efficiency of ~0.5. These movies should
%   include the full camera sensor.
%
%   calcLocalGamma(...) without an output arguments will prompt the user
%   to save the mapping as a .mat file.
%
%   calcLocalGamma() without input arguments will prompt for files.
% 
% See also: applyLocalGamma, calc_gamma, gammacorrect, correctTraces.

% Copyright 2025 St Jude Children's Research Hospital. All Rights Reserved.

% FIXME:
%  - load input data only once, keeping only mean fret values and coordinates.


% PARAMETERS
maskCorners   = true;  %ignore field of view corners with more aberration.
contourStride = 30;    %pixel stride length of FRET variance contour plot
contourFrames = 1:25;  %number of frames to average to get FRET value
minFret = 0.3;         %reject traces with data outside this range.
maxFret = 0.7;         %...adjust for different standard samples.
contourLevels = 0.7:0.05:1.3;  %levels for display of contour plots
fretaxis = 0:0.02:1;   %FRET histograms levels

warning('OFF','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');


% Choose gamma-corrected autotrace files
if nargin<1
    result = [];
    files = getFiles('*.traces','Select calibration data',false);
    if numel(files)==0, return; end  %user hit cancel
elseif ischar(files) || isstring(files)
    files = {files};
else
    assert(iscell(files),'Invalid input. Should be path to traces file or cell array of paths.');
end
nFiles = numel(files);


% Calculate mean FRET ratio (A/D) over field of view from all input files.
[X,Y,Z] = deal(zeros(0,1));
[nX,nY] = deal(0);  %field of view size in pixels

for i=1:nFiles
    data = loadTraces(files{i});
    fret = data.fret(:,contourFrames);
    selected = all( fret>minFret & fret<maxFret, 2);
    data.subset(selected);
    fprintf('Selected %d/%d (%0.1f%%) traces\n',...
            sum(selected), numel(selected), 100*sum(selected)/numel(selected));

    fret = data.acceptor(:,contourFrames) ./ data.donor(:,contourFrames);
    Z = [Z double(mean(fret,2)/mean(fret(:)))']; %#ok<*AGROW> 
    X = [X [data.traceMetadata.donor_x]];
    Y = [Y [data.traceMetadata.donor_y]];

    if all(isfield(data.fileMetadata,{'nX','nY'}))
        nX = max(nX,data.fileMetadata.nX);
        nY = max(nY,data.fileMetadata.nY);
    end
end


% Interpolate values to create a contour plot of field of view
% Contour plot of obserevd mean FRET over field of view
figure;
subplot(1,4,1);

[Xq,Yq] = meshgrid(0:contourStride:nX, 0:contourStride:nY);
Zq = griddata(X,Y,Z, Xq,Yq, 'natural');
contourf(Xq,Yq,Zq, contourLevels, 'LineColor','none');
clim( contourLevels([1 end]) );
colorbar;
axis square;
ylabel('Y Position (px)');
xlabel('X Position (px)');
title('Data');


% Fit above to a 2nd order polynomial and display as contour plot.
% Corners are masked, which might not be relevant to all datasets!
if maskCorners
    roi = images.roi.Circle('Center',[nX/2,nY/2], 'Radius',max(nX,nY)/2);
    selected = roi.inROI(X,Y);
    X = X(selected);
    Y = Y(selected);
    Z = Z(selected);
end
result = fit([X' Y'], Z', 'poly22');

subplot(1,4,2);
contourf(Xq,Yq,result(Xq,Yq), contourLevels, 'LineColor','none');
clim( contourLevels([1 end]) );
c = colorbar;
c.Label.String = '<A/D>';
axis square;
xlabel('X Position (px)');
title('Model')


% Correct the calibration data as a sanity check (should get sharper)
[rawhist,corrhist] = deal( zeros(1,numel(fretaxis)-1) );
[X,Y,Z] = deal(zeros(0,1));

for i=1:nFiles
    data = loadTraces(files{i});

    % Mask out the corners (optional)
    if maskCorners
        x = [data.traceMetadata.donor_x]';
        y = [data.traceMetadata.donor_y]';
        roi = images.roi.Circle('Center',[nX/2,nY/2], 'Radius',max(nX,nY)/2);
        data.subset( roi.inROI(x,y) );
    end

    % Histograms before and after correction
    rawhist = rawhist + histcounts( data.fret(:,contourFrames), fretaxis );

    x = [data.traceMetadata.donor_x]';
    y = [data.traceMetadata.donor_y]';
    data.acceptor = data.acceptor ./ result(x,y);
    data.recalculateFret();
    corrhist = corrhist + histcounts( data.fret(:,contourFrames), fretaxis );

    % Residual "error"
    fret = data.fret(:,contourFrames);
    Z = [Z double(mean(fret,2)/mean(fret(:)))'];
    X = [X [data.traceMetadata.donor_x]];
    Y = [Y [data.traceMetadata.donor_y]];
end

rawhist = 100*rawhist/sum(rawhist);
corrhist = 100*corrhist/sum(corrhist);


% Corrected FRET contour plot to evaluate residual "error".
subplot(1,4,3);
Zq = griddata(X,Y,Z, Xq,Yq, 'natural');
contourf(Xq,Yq,Zq, contourLevels, 'LineColor','none');
clim( contourLevels([1 end]) );
colorbar;
axis square;
ylabel('Y Position (px)');
xlabel('X Position (px)');
title('Residual');
axis square;

% Fit distributions to Gaussians for a quantitative comparison.
rawModel = fit( fretaxis(1:end-1)', rawhist', 'Gauss1' );
ffModel  = fit( fretaxis(1:end-1)', corrhist', 'Gauss1' );
fprintf('Input  = %.3f +/- %.3f\n', rawModel.b1, rawModel.c1/sqrt(2));
fprintf('Output = %.3f +/- %.3f\n', ffModel.b1, ffModel.c1/sqrt(2));

% Overlay FRET histograms before and after correction.
subplot(1,4,4);
plot(fretaxis(1:end-1),rawhist','ob');
hold on;
plot(fretaxis(1:end-1),corrhist','or');
plot(rawModel,'b-');
plot(ffModel,'r-');
legend({'raw','corrected'});
xlim([0.3 0.7]);
xlabel('FRET');
ylabel('Counts (%)');
axis square;


if nargout==0
    [f,p] = uiputfile('*.mat','Save correction model','flatfield.mat');
    if isequal(f,0), return; end
    save(fullfile(p,f),'nX','nY','result');
end


end %function

