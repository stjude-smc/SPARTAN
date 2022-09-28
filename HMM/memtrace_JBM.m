function memtrace_JBM(input)
%
%

narginchk(0,2);

MIN_DWELLS = 20;  % minimum number of dwells per trace


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


% Plot autocorrelation for each state.
temp = cat(1,dwt{:});
nStates = max( temp(:,1) );
figure;

for i=1:nStates
    subplot(1,nStates,i);
    dwtcorr(dwt, i, MIN_DWELLS);
    title( sprintf('State %d',i) );
end


end




%% Calculate the autocovariance of the dwell times in each trace
function dwtcorr(dwt, STATE, MIN_DWELLS)

ntraces = numel(dwt);
sel = false(ntraces,1);

lens_ac = zeros(ntraces,1);
autc = cell(ntraces,1);

% Calculate dwell-time autocorrelation for traces w/ minimum num. dwells
for ntr = 1:ntraces 
    dwells = dwt{ntr};
    dwells = dwells( dwells(:,1)==STATE, : );
    
    if size(dwells,1)<MIN_DWELLS, continue; end
    sel(ntr) = true;
    
    % Calculate the autocovariance in the dwell time sequence
    ac = acorrdwt(dwells);
    autc{ntr} = ac;
    lens_ac(ntr) = numel(ac);
end

% Calculate average autocovariance at each time point
lens1 = lens_ac(sel);
tmp1 = autc(sel,1);
a1 = max(lens1);
ac1 = zeros(numel(tmp1),a1);

% Make matrices with traces in rows, dwell indeces in columns
for i = 1:numel(tmp1)
    ac1(i,1:lens1(i)) = tmp1{i,1}';  
end
    
avgac1 = zeros(1,a1);
acerrbars1 = zeros(1,a1);

for i = 1:a1
    n = sum(lens1 >= i); %counts how many ac arrays have i+ elements (i+ dwells)
    if n <= 1            %if there are 2+ arrays with length i+, do stats on them
       break 
    end
    avgac1(i) = nanmean(ac1(lens1 >= i,i)); %takes mean & sem over column i (dwell index i) and all rows with i+ length (lens1 >= i)
    acerrbars1(i) = nanstd(ac1(lens1 >= i,i))/sqrt(n);
end

temp_s1 = fitoptions('Method','NonlinearLeastSquares','Lower',[-0.1 0 0],'Upper',[0.01 1 1000],...
        'Startpoint',[0 0.06 5]);
temp_f1 = fittype('y0+A1*exp(-k1*(x-1))','dependent',{'y'},'independent',{'x'},...
        'coefficients',{'y0','A1','k1'}, 'options',temp_s1);

time1 = 1:numel(avgac1);
result1=fit(time1(2:50)', avgac1(2:50)', temp_f1); 

plot(time1(1:end),avgac1(1:end),'ok','MarkerFaceColor','b')
hold on
errorbar(time1(1:end),avgac1(1:end),acerrbars1(1:end),'ok','MarkerFaceColor','b')
plot( result1(1:50), 'r-')
xlabel('Dwell number');
ylabel('Autocorrelation');
grid on
zoom on

ap = get(gca,'Position');  %left bottom width height
annotation( gcf, 'textbox', [ap(1)+0.9*ap(3) ap(2)+0.9*ap(4) 0.1*ap(3) 0.1*ap(4)], ...
            'String',sprintf('N=%d', sum(sel)), 'HorizontalAlignment','right', ...
            'LineStyle','none', 'tag','Nmol' );

% axis([0 50 -0.1 0.25])  %for simulations
axis([0 50 -0.05 0.1])

end



function ac = acorrdwt(dwt)
% Calculates the autocorrelation of a sequence of dwells times.
% Uses Matlab xcov function (Signal Processing Toolbox)

dwells = dwt(:,2);
ndw = numel(dwells);

temp = xcov(dwells,'coeff');
ac = temp(ndw:end);

end
