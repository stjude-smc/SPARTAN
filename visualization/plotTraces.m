function plotTraces(filename,plotTitle)
% PLOTTRACES  Plots fluorescence and idealized FRET traces
%
%   plotTraces() ...

%   Copyright 2007-2015 Cornell University All Rights Reserved.


gamma = 1;
allowMinutes = true;

% If not files specified, prompt user for them.
if ~exist('filename','var'),
    [datafile,datapath] = uigetfile({'*.traces;*.txt'},'Choose a traces file:');
    if datafile==0, return; end  %user hit "cancel"
    filename = [datapath filesep datafile];
end


% Generate plot titles 
if ~exist('plotTitle','var'),
    [p,name] = fileparts(filename);
    plotTitle = strrep(name,'_',' ');
end


% Load fluorescence and FRET data
% [donor,acceptor,fret,ids,time] = loadTraces( filename );
data = loadTraces( filename );
data.fret(data.fret>1) = 1;
data.fret(data.fret<-0.5) = -0.5;
time = data.time;
[nTraces,len] = size(data.donor);

inFrames = (time(1)==1);
if ~inFrames
    sampling = time(2)-time(1)
    time = time/1000; %in seconds
else
    f = inputdlg('What is the sampling interval (in ms) for this data?');
    sampling = str2double(f)
    time = (sampling/1000).*(0:len-1);
end

% Load idealization data
[p,n] = fileparts(filename);
dwt_fname = fullfile(p, [n '.qub.dwt']);
if exist(dwt_fname,'file')
    [dwt,dwtSampling,offsets,model] = loadDWT( dwt_fname );
    model = model(:,1);
    dwtSampling = double(dwtSampling)
    
    % Look for simulated idealization data
    if nTraces==length(dwt),
        dwtToTraceID = 1:nTraces;
    else
        dwtToTraceID = floor( (offsets+1)/len )+1;
    end
    
    if sampling ~= dwtSampling,
        warning('plotTraces: data and idealization sampling intervals do not match');
    end
else
    dwtToTraceID = [];
end



%% Create GUI environment for cycling through rates


f1 = figure;

set(f1,'defaultaxesfontsize',16);
set(f1,'defaulttextfontsize',18);
set(f1,'Position',[264   405   672   513]);


%left bottom width height
slider = uicontrol( f1, 'Style','slider', 'Callback',@showByIndex, ...
    'Min',1, 'Max',nTraces, 'SliderStep',[1/nTraces 10/nTraces], 'Value',1, ...
    'Units','normalized', 'Position',[.933 .33 .04 .4] );

molecule_no = 1;
showTrace;  %display initial view
    



function showPrev(h, eventdata)
    molecule_no = molecule_no-1;
    if molecule_no==1,
        set(prev, 'Enable', 'off');
    elseif molecule_no<nTraces
        set(next, 'Enable', 'on');
    end

    showTrace();
end  %function ShowPrev

function showNext(h, eventdata)
    molecule_no = molecule_no+1;
    if molecule_no==nTraces,
        set(next, 'Enable', 'off');
    elseif molecule_no>1
        set(prev, 'Enable', 'on');
    end

    showTrace();
end  %function ShowNext

function showByIndex(h, eventdata)
    molecule_no = floor( get(slider,'Value') );
    showTrace();
    
    get(gca,'Position');
end

function showTrace

    if allowMinutes && sampling>100,
        sampFact=60;
    else
        sampFact=1;
    end

    %constants
    nrows = 2;
    ncols = 1;
    nBins = 60;
    
    ax = [];
    
    %
    i = molecule_no;
    dwtID = find( dwtToTraceID==i, 1 );
    

    % Load convert DWT to idealization
    if ~isempty(dwtID),
        states = dwt{dwtID}(:,1);
        times  = double(dwt{dwtID}(:,2));
        dwtTotalTime = sum(times);
        sidl = [];
        
        for j=1:numel(states),
            dtime = times(j);
            sidl = [sidl repmat( states(j), 1,dtime ) ];
        end
        
        newModel = model;
%         for j=1:max(states),
%             newModel(j) = mean( fret(i,sidl==j) );
%         end
        
        %disp(newModel)
        idl = newModel( sidl );
        
        %idl = [ idl repmat(model(1),1,len-dwtTotalTime) ];
        
        dwtTime = (dwtSampling/1000).*(0:dwtTotalTime-1);
    end
    % Draw fluorescence time trace
    ax(1) = subplot(nrows,ncols,1);

    d = gamma*data.donor(i,:) /1000;
    a = data.acceptor(i,:) /1000;
    
%     sd  = sort(d);  sa = sort(a);
%     sdm = sd( end-ceil(0.01*length(sd)) );
%     sam = sa( end-ceil(0.01*length(sa)) );
    
    plot( time/sampFact, d,'g', 'LineWidth',1  ); hold on;
    plot( time/sampFact, a,'r', 'LineWidth',1  ); hold on;
    ylabel('Fluorescence');
    ymax = max( max(d)+3,max(a)+3 );
    ylim( [-3 ymax] );
    set(gca,'xticklabel',[]);
    grid on;
    zoom on;
    hold off;
    title( [plotTitle ' - ' num2str(molecule_no)] );
%     xlabel( 'Time (sec)' );

    % Draw FRET time trace
    ax(2) = subplot(nrows,ncols,ncols+1);
    cla;
    
    plot( time/sampFact, data.fret(i,:), 'b-', 'LineWidth',1 );

    % Draw idealization on top of FRET trace
    if ~isempty(dwtID),
        hold on;
        stairs( dwtTime/sampFact, idl, 'r-', 'LineWidth',1 );
    end
    
    grid on;
    zoom on;
    ylabel('FRET');
    ylim( [-0.1 1.0] );


    pbTime = (find( data.fret(i,:)~=0, 1, 'last' )*1.04) * sampling/sampFact/1000;

%     maxlife = sum(times) * sampling/1000 +0.5;
    maxlife = max( pbTime, 5 );
    xlim( [0 maxlife] );
    if allowMinutes && sampling>100
        xlabel( 'Time (min)' );
    else
        xlabel( 'Time (sec)' );
    end
    set(gca,'ytick', [0, 0.2, 0.4, 0.6, 0.8, 1.0]);
    hold off;

    
    %subplot(2,1,1);  set(gca,'xticklabel',[]);
    
    
    %--- Plot distributions of FRET and fluorescence in the highlighted window
    
    if ncols>=2,
    
        % Generate the histograms
        pbFrame = find( data.fret(i,:)~=0, 1, 'last' )-1;

        f = [data.acceptor(i,:) data.donor(i,:)]/1000;
        fluor_bins = min(f)-2:(range(f)/nBins):max(f)+1;

        f = data.fret(i,:);
        fret_bins = min(f)-0.1:(range(f)/nBins):max(f)+0.1;

        [ha] = hist( data.acceptor(i,1:pbFrame)/1000, fluor_bins);
        [hd] = hist( data.donor(i,1:pbFrame)/1000, fluor_bins);
        [hf] = hist( data.fret(i,1:pbFrame), fret_bins);

        % Plot distribution of fluorescence
        ax(3) = subplot(nrows,ncols,2);

        plot( ha, fluor_bins, 'r-', 'LineWidth',2 );
        hold on;
        plot( hd, fluor_bins, 'g-', 'LineWidth',2 );
        hold off;
        grid on;
        set(gca,'xticklabel',[]);
        pbaspect( [1 1 1 ] );

        yax = [ax(1) ax(3)];
        linkaxes(yax,'y');

        % Plot distribution of FRET
        ax(4) = subplot(nrows,ncols,4);

        stairs( hf, fret_bins, 'b-', 'LineWidth',2 );
        pbaspect( [1 1 1] );
        grid on;
    %     set(ax(4),'yticklabels',[]);

        yax = [ax(2) ax(4)];
        linkaxes(yax,'y');
    end
    
        
    % Hold the time axis consistent
    linkaxes([ax(1) ax(2)],'x');
    
    xlim( ax(1), [0 maxlife] );
    
%     xlim([0 10])

    if nrows>=3,
        ax3 = subplot(3,1,3);
        
        mfd = medianfilter(data.donor(i,:),9);
        mfa = medianfilter(data.acceptor(i,:),9);
        
        s = -gradient(mfd).*gradient(mfa);
        s(s<0)=0;
        
        plot( time/sampFact, s );
        grid on;
          
        linkaxes([ax(1) ax(2) ax3],'x');
    end
    

end %function showTrace


end %function plotTraces
