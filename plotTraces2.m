function plotTraces2(filename,plotTitle)
% PLOTTRACES  Plots fluorescence and idealized FRET traces
%
%   plotTraces() ...
%


gamma = 1;
% model = [0.01 0.25 0.39 0.58];

% If not files specified, prompt user for them.
if ~exist('filename','var'),
    [datafile,datapath] = uigetfile({'*.txt'},'Choose a traces file:');
    if datafile==0, return; end  %user hit "cancel"
    filename = [datapath filesep datafile];
end


% Generate plot titles 
if ~exist('plotTitle','var'),
    [p,name] = fileparts(filename);
    plotTitle = strrep(name,'_',' ');
end


% Load fluorescence and FRET data
[donor,acceptor,fret] = LoadTraces( filename );
[nTraces,len] = size(donor);

% Load idealization data
dwt_fname = strrep( filename, '.txt', '.qub.dwt' );
if exist(dwt_fname,'file')
    [dwt,framerate,offsets,model] = LoadDWT( dwt_fname );
    model = model(1:2:end);
    framerate = double(framerate);
    dwtToTraceID = floor( (offsets+1)/len )+1;
end




%% Create GUI environment for cycling through rates


f1 = figure;

set(f1,'defaultaxesfontsize',16);
set(f1,'defaulttextfontsize',18);
set(f1,'Position',[204 331 672 577]);

%left bottom width height
slider = uicontrol( f1, 'Style','slider', 'Callback',@showByIndex, ...
    'Min',1, 'Max',nTraces, 'SliderStep',[1/nTraces 10/nTraces], 'Value',1, ...
    'Units','normalized', 'Position',[.933 .33 .04 .4] );

% prev = uicontrol(f1, 'String', 'Prev', 'Callback', @showPrev, ...
%     'Units', 'pixels', 'Position', [20 10 80 30]);
% next = uicontrol(f1, 'String', 'Next', 'Callback', @showNext, ...
%     'Units', 'pixels', 'Position', [100 10 80 30]);

molecule_no = 1;
% set(prev, 'Enable', 'off');
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

    i = molecule_no;
    dwtID = find( dwtToTraceID==i, 1 );
    
    time = (framerate/1000).*(0:len-1);

    % Load convert DWT to idealization
    if ~isempty(dwtID),
        states = dwt{dwtID}(:,1);
        times  = double(dwt{dwtID}(:,2));
        tlen = sum(times);
        idl = [];

        for j=1:numel(states),
            dtime = times(j);
            idl = [idl repmat( model(states(j)), 1,dtime ) ];
        end
        idl = [ idl repmat(model(1),1,len-tlen) ];
    end

    % Draw fluorescence time trace
    ax(1) = subplot(2,1,1);

    d = gamma*donor(i,:) /1000;
    a = acceptor(i,:) /1000;
    plot( time, d,'g', time,a,'r' );
    ylabel('Fluorescence');
    ymax = max( max(a)+3,max(d)+3 );
    ylim( [-3 ymax] );
    set(gca,'xticklabel',[]);
    grid on;
    zoom on;
    hold off;
    title( [plotTitle ' - ' num2str(molecule_no)] );

    % Draw FRET time trace
    ax(2) = subplot(2,1,2);

    plot( time, fret(i,:), 'b' );
    ylabel('FRET');
    grid on;
    zoom on;
    ylim( [-0.1 1.0] );

    % Draw idealization on top of FRET trace
    if ~isempty(dwtID),
        hold on;
        stairs( time, idl, 'r-' );
    end

    pbTime = (find( fret(i,:)~=0, 1, 'last' )+10) * framerate/1000;
%     maxlife = sum(times) * framerate/1000 +0.5;
    maxlife = max( pbTime, 5 );
    xlim( [0 maxlife] );
    xlabel( 'Time (s)' );
    set(gca,'ytick', [0, 0.2, 0.4, 0.6, 0.8, 1.0]);
    hold off;
    
    % Hold the time axis consistent
    linkaxes(ax,'x');

end %function showTrace




end %function plotTraces
