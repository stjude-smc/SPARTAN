function [output,binCenters,labels] = ratehist( rates, edges )
% ratehist: display histograms of rate constants from an array of models.
%
%   HIST = ratehist( RATES ) calculates histograms as columns
%   in HIST of the rate constants in RATES, which is a 3-dimensional matrix
%   with the rate matrix in the first two dimensions and traces across the
%   third dimension. The first column of HIST contains histogram bin edges.
%   The histogram is on a log scale.
%
%   ratehist( RATES ) without output arguments plots the
%   histogram data in a new figure instead.


% Check input parameters
narginchk(1,2);
nargoutchk(0,3);

assert( isnumeric(rates), 'Input must be a matrix of rates' );
if nargin<2
    edges = -5:0.2:5;
end


%% Calculate histograms

% Set order of state pairs that describe each rate constant
[src,dst] = find( all(rates>0,3) );  %& ~model.fixRates;
[src,idx] = sort(src);
dst = dst(idx);

binCenters = ( edges(1:end-1) + edges(2:end) )/2;
output = zeros( numel(binCenters), numel(src) );
labels = zeros( 0, 2 );

% Lump extreme rates into first and last bins.
% Otherwise they will not be represented in the histograms!
logrates = log10( rates );
logrates( logrates>=edges(end) ) = edges(end)-0.01;
logrates( logrates<=edges(1) ) = edges(1)+0.01;

% Create a log-scale histogram for each rate constant in the model
if nargout<1
    hFig = figure;
    ax = zeros( numel(src), 1 );
end

for i=1:numel(src)
    % Calculate histogram for this rate (i->j)
    values = logrates( src(i), dst(i), : );
    values = values(~isnan(values));
    output(:,i) = histcounts( values, edges );
    labels(i,:) = [src(i), dst(i)];

    % Display the histogram plots (if requested)
    if nargout<1
        ax(i) = subplot( 1, numel(src), i );
        bar( ax(i), binCenters, output(:,i) );
        title( ax(i), sprintf('%d -> %d',src(i),dst(i)) );
        xticks( ax(i), floor(edges(1)):1:ceil(edges(end)) );
        xlabel('log10(rate)');
        if i==1, ylabel('Traces'); end
    end
end

if nargout<1
    linkaxes( ax, 'x' );
    
%     %% Add menu items for adjusting settings and saving output to file
%     prompt = {'Remove blinks:', 'Log scale:', 'Log bin size:', 'Normalization:'};
%     fields = {'removeBlinks', 'logX', 'dx', 'normalize'};
%     types{4} = {'none','state','file','time'};
%     output = [to_col(dwellaxis) horzcat(histograms{:})];
%     outFit = [to_col(fitaxis) fits];
% 
%     defaultFigLayout( hFig, @(~,~)ratehist(getFiles('*.dwt'),params),    ...  %File->New
%                             @(~,~)ratehist(hFig,getFiles('*.dwt'),params), ...  %File->Open
%                            {@exportTxt,dwtfilename,output}, ...
%            {'Copy histograms',{@clipboardmat,output} );
    defaultFigLayout( hFig, [],    ...  %File->New
                            [], ...  %File->Open
                            [], ...  %Export
           {'Copy histograms',{@clipboardmat,[binCenters',output]}} );  %Edit
end





end  %function



