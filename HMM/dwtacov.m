function varargout = dwtacov(input)
% dwtacov  autocovariance of dwell times in each state, mean over traces
%
%   dwtacorr(INPUT) plots the autocovariance of dwell times in each state 
%   and fits to a single exponential. INPUT can be:
%   - matrix: state assignment at each frame (columns) and trace (rows).
%   - cell array of traces: matrix with state (1st col) and time (2nd col).
%   - string: path to a .dwt file.
%
%   dwtacorr() with no input arguments will prompt for a .dwt file path.
%
%   [AC,ACERR] = dwtacorr(...) returns the autocorrelation curves (AC) with
%   states in columns and lag times over rows and standard errors of these
%   values (ACERR) over the list of traces.
%
%   Requires Signal Processing Toolkit (for xcov).

% Copyright 2022, Cornell University and St Jude Childrens Research Hospital

narginchk(0,1);

% minimum number of dwells per trace reqquired to calculate autocorrelation. 
% without this, early lags may have far more traces than later lags, which
% can significantly distort the curve.
MIN_DWELLS = 20;

% Process input arguments to get dwell-time array
if nargin<1
    [f,p] = uigetfile('*.dwt','Choose dwt file:');
    if f==0, return, end
    dwt = loadDWT( fullfile(p,f) );
elseif ischar(input)
    dwt = loadDWT(input);
elseif isnumeric(input)
    dwt = idlToDwt(input);
else
    error('Invalid first input');
end

temp = cat(1,dwt{:});
nStates = max( temp(:,1) );
[ac,acerr] = deal( cell(nStates,1) );

% Plot autocorrelation for each state.
figure;
ax = zeros(nStates,1);
for i=1:nStates
    ax(i) = subplot(1,nStates,i);
    [ac{i},acerr{i}] = calc_dwtacorr(dwt, i, MIN_DWELLS);
    title( sprintf('State %d',i) );
end
linkaxes(ax,'xy');

% Prepare output arguments
ac = [ac{:}];
acerr = [acerr{:}];
output = {ac,acerr};
[varargout{1:nargout}] = output{1:nargout};

end



%% Calculate the autocovariance of the dwell times in each trace
function [avgac,acerrbars] = calc_dwtacorr(dwt, STATE, MIN_DWELLS)

maxlag = 50;
ntraces = numel(dwt);
sel = false(ntraces,1);
ac1 = nan(ntraces,maxlag);

% Calculate autocovariance curve for each trace
for ntr = 1:ntraces 
    dwells = dwt{ntr};
    dwells = dwells( dwells(:,1)==STATE, 2 );
    ndwell = numel(dwells);

    % If there aren't MIN_DWELLS dwells, skip the trace.
    if ndwell<=MIN_DWELLS, continue; end
    sel(ntr) = true;
    
    % Calculate the autocovariance in the dwell time sequence
    %[ac,lag] = xcorr( dwells-mean(dwells), maxlag, 'coeff' );
    [ac,lag] = xcov(dwells,maxlag,'coeff');
    ac = ac( lag>=0 & lag<=ndwell );
    ac1(ntr,1:numel(ac)) = ac; 
end

% Average covariance curves over traces
avgac = zeros(maxlag,1);
acerrbars = zeros(maxlag,1);

for i = 1:maxlag
    % Covariance at lag i from traces with at least that many dwells.
    data = ac1( sel, i );
    data = data( ~isnan(data) );
    
    if numel(data)<3, break; end
    
    avgac(i) = mean(data);
    acerrbars(i) = std(data)/sqrt(numel(data));
end

% Fit autocorrelation to an exponential function.
temp_s1 = fitoptions('Method','NonlinearLeastSquares','Lower',[-0.1 0 0],'Upper',[0.01 1 1000],...
        'Startpoint',[0 0.06 5]);
temp_f1 = fittype('y0+A1*exp(-k1*(x-1))','dependent',{'y'},'independent',{'x'},...
        'coefficients',{'y0','A1','k1'}, 'options',temp_s1);
time = 1:numel(avgac);
result = fit(time(2:end)', avgac(2:end), temp_f1); 

% Overlay autocorrelation and fit
plot(time,avgac,'ok','MarkerFaceColor','b');
hold on
errorbar(time,avgac,acerrbars,'ok','MarkerFaceColor','b');
plot( result(1:maxlag), 'r-');
xlabel('Dwell number');
ylabel('Autocorrelation');
grid on;
zoom on;

ap = get(gca,'Position');  %left bottom width height
annotation( gcf, 'textbox', [ap(1)+0.9*ap(3) ap(2)+0.9*ap(4) 0.1*ap(3) 0.1*ap(4)], ...
            'String',sprintf('N=%d', sum(sel)), 'HorizontalAlignment','right', ...
            'LineStyle','none', 'tag','Nmol' );

% axis([0 50 -0.1 0.25])  %for simulations
axis([0 50 -0.05 0.1])


% Add menus to change settings, get data, open new plots, etc.
% fields = {'MIN_DWELLS',             'maxlag'};
% prompt = {'Min. number of dwells:', 'Max. lag time'};
% 
% defaultFigLayout( hFig, @(~,~)frethistComparison(getFiles(),params), ...  %File->New callback
%                         @(~,~)frethistComparison(cax,getFiles(),params), ...  %File->Open callaback
%                         {@exportTxt,files,output}, ...                     %Export callaback
%       { 'Change settings...',@(~,~)settingdlg(params,fields,prompt,@frethistComparison,{cax,files}); ...
%         'Reset settings',    @(~,~)frethistComparison(cax,files); ...
%         'Copy values',      {@clipboardmat,output}  }  );


end

