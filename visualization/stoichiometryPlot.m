function varargout = stoichiometryPlot(varargin)
% stoichiometryPlot  FRET vs. stoichiometry contour plot (all frames)
%
%   stoichiometryPlot( data )
%   Draws the transition density contour plot in TDP.  If no arguments
%   supplied, the function will prompt the user for the file.
%
%   stoichiometryPlot( AX, ... )
%   Draws the plot in the axes object AX.

%   Copyright 2022 All Rights Reserved.


% Extract axes if specified as the first argument. "args" are the remaining
% arguments.
[cax,args] = axescheck(varargin{:});
isNewAx = isempty(cax);

if numel(args)<1,
    error('Not enough input arguments');
end


% Load optional arguments
options = cascadeConstants('defaultMakeplotsOptions');

if numel(args)==2,
    assert( isstruct(args{2}) );
    options = mergestruct( options, args{2} );

elseif numel(args)>2,
    args = args(2:end);
    assert( iscell(args) & mod(numel(args),2)==0, ...
            'Incorrect format for optional arguments list' );
    vopt = struct(args{:});
    options = mergestruct( options, vopt );
end

    
% Load trace data
if ischar(args{1}) || isstring(args{1})
    data = loadTraces(args{1});
elseif isa(args{1},'Traces')
    data = args{1};
else
    error('Invalid input argument: should be a Traces object or path to a .traces file');
end
assert( data.isChannel('acceptorDirect'), 'This is not ALEX data!' );


% Create histogram
frames = options.pophist_offset+(1:options.contour_length);
st = data.stoichiometry(:,frames);
fret = data.(options.fretField)(:,frames);

% Remove periods where donor is dark/bleached.
st = st(fret~=0);
fret = fret(fret~=0);

counts = histcounts2( fret, st, options.fret_axis, options.fret_axis );
counts = counts / sum(counts(:));

% FIXME: may need to pad with NaN if there is just one trace??



%% Draw the filled contour plot in current axis
max_mol = 1/(5*options.cplot_scale_factor);
nl = size(options.cmap,1)-1;       %number of contour levels
con = 0:(max_mol/nl):max_mol;      %contour levels

if isNewAx, figure; end
cax = newplot(cax);

binCenters = options.fret_axis(1:end-1)+diff(options.fret_axis(1:2))/2;
[~,hPlot] = contourf(cax, binCenters, binCenters, counts', con);
colormap(cax, options.cmap);
set(hPlot, 'LineColor', 'none');

% Add all appearance details
xlabel(cax,'FRET');
ylabel(cax,'Stoichiometry');
zoom(cax,'on');
set(cax,'xtick', 0:0.2:1, 'ytick', 0:0.2:1);
set(cax,'xlim',[-0.1 1.1]);
set(cax,'ylim',[-0.1 1.1]);
% axis(cax, bounds);
set(cax, 'XGrid','on', 'YGrid','on', 'Box','on', 'PlotBoxAspectRatio',[1 1 1]);

output = {counts'};
[varargout{1:nargout}] = output{:};



end

