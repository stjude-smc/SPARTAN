function tplot(varargin)
% TPLOT  Draws a transition density contour plot
%
%   TPLOT( TDP, ... )
%   Draws the transition density contour plot in TDP.  If no arguments
%   supplied, the function will prompt the user for the file.

%   Copyright 2007-2015 Cornell University All Rights Reserved.

%Plots transition density plot created by tdplot from a file, or from the
%variable tdp entered at command line.


% Extract axes if specified as the first argument. "args" are the remaining
% arguments.
[cax,args] = axescheck(varargin{:});
cax = newplot(cax);

% Load TD Plot data
if numel(args)<1,
    error('Not enough input arguments');
end

tdp = args{1};
assert( isnumeric(tdp) );


% Load optional arguments
constants = cascadeConstants();
options = constants.defaultMakeplotsOptions;

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


f_axis=tdp(2:end,1);
i_axis=tdp(1,2:end);

% Setup contour levels
top = options.tdp_max;
con=0:(top/size(options.cmap,1)):top;
tdpfix = tdp;
tdpfix(end,end) = 10;  % hack to make colorscale fixed

% draw contour plot
[~,hand]=contourf(cax,i_axis,f_axis,tdpfix(2:end,2:end),con);

% Extra formatting
set(hand,'LineColor','none');
colormap(cax,options.cmap);
set(cax, 'PlotBoxAspectRatio', [1 1 1]);
ylim(cax, options.fretRange );
xlim(cax, options.fretRange );
ylabel(cax,'Final FRET');
xlabel(cax,'Inital FRET');
set(cax, 'XGrid','on', 'YGrid','on', 'Box','on');


% Display number of transitions and transition rate.
if ~options.hideText,
    total_time = tdp(1,1);
    t = sum(sum( tdp(2:end,2:end)*total_time ));

    textOpt = {'HorizontalAlignment','center', 'Parent',cax};
    text( 0.45,0.9, sprintf('N_t=%.0f',t),            textOpt{:} );
    text( 0.45,0.0, sprintf('t/s=%.2f',t/total_time), textOpt{:} );
end

