%% MAKEFIGURE: constructs statehist histrograms for figure preparation
function makefigure( filenames )

constants = cascadeConstants();
darkModel = constants.modelLocation;
sampling = 40; %doesn't really matter, as long as its close
contour_bin_size = 0.02;
sumlen = 50;
fretaxis = -0.1:contour_bin_size:1.0;
nbins = length(fretaxis);


% Select data files to display from user
if nargin<1,
    filenames = getFiles;
    if numel(filenames)==0, return; end
end

nFiles = numel(filenames);


% Load each data file and idealize to a 2-state model
% and build histograms using only the non-zero FRET state.
frethist = zeros(nbins,nFiles);

model = qub_loadModel( [darkModel 'darkstate.qmf'] );
model.fixMu    = ones(1,2);
model.fixSigma = ones(1,2);
model.mu(2) = 0.3;
skmParams.quiet = 1;

for i=1:nFiles
    
    % Load FRET data from file
    [d,a,data] = loadTraces( filenames{i} );
    [nTraces,nFrames] = size(data);
    
    % Idealize to two-state model to select dark state
    [dwt,newModel] = skm( data, sampling, model, skmParams );
    
    % Expand DWT into an idealization: 2=non-zero FRET state
    idl = zeros(nTraces,nFrames);
    
    for dwtID=1:nTraces,
        states = dwt{dwtID}(:,1);
        times  = double(dwt{dwtID}(:,2));

        ends = cumsum( times );
        starts = cumsum( [1; times(1:end-1)] );
        for j=1:numel(states),
            idl(dwtID,starts(j):ends(j)) = states(j);
        end
    end
    
    % Find all datapoints in non-zero FRET state and
    % Build a histogram of the data and save it
    fh = hist( data(idl==2), fretaxis )';
    frethist(:,i) = fh/sum(fh);
    
end


%s=get(0,'ScreenSize');
%figure('Position',s);
f1 = figure;
set(f1,'defaultaxesfontsize',16);
set(f1,'defaulttextfontsize',18);

colors=[0   0   0;...
        1   0   0;...
        1   0   0;...
        1   0   0];

cla;
fillStairs( fretaxis, frethist, colors );

% [file,path]=uiputfile('*.jpg','Save figure image as:');
% outfile=strcat(path,file);
% imwrite(gcf,outfile,'jpg','Quality',100)


end







function fillStairs( x, sthist, colors )

binsize=x(2)-x(1);


% % Draw solid outline of total histogram
% % stairs(x(26:end)-bin/2,100*sum(sthist(26:end,2:end),2),'Color',[0.7 0.7 0.7])
% totalhist = sum(sthist,2)*100;
% stairs(x-binsize/2,totalhist,'Color',[0.7 0.7 0.7]);
% hold on;
    
    
[nBins,nSamples] = size(sthist);
    
for j=1:nSamples

    y=sthist(:,j)*100;

    xx=zeros(2*length(x),1);
    yy=zeros(size(xx));
    for i=1:2:2*length(x)
        xx(i)=x((i+1)/2);
        xx(i+1)=x((i+1)/2)+binsize;
        yy(i)=y((i+1)/2);
        yy(i+1)=y((i+1)/2);
    end
    
    xx = [xx; xx(end)+binsize];
    yy = [yy; 0];

    if j==3 || j==4
        patch(xx-binsize/2,yy,colors(j,:),'EdgeColor','none','FaceAlpha',0.25);
    else
        patch(xx-binsize/2,yy,colors(j,:),'EdgeColor','none','FaceAlpha',0.15);
    end
    hold on;
    stairs(x-binsize/2,y,'Color',colors(j,:));
    hold on;

end


axis([-0.05 1.001 0 8])
% set(gca,'PlotBoxAspectRatio',[1.7 1 1]);
ylabel('Probability (%)')
xlabel('FRET')
hold off

end















