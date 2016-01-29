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
    [tdpfile, tdppath]=uigetfile('*_tdp.txt','Choose tdplot:');
    if tdpfile==0
        disp('No File Selected.')
        return;
    else
        file=strcat(tdppath,tdpfile);
        tdp=load(file);
    end
else
    tdp = args{1};
    assert( isnumeric(tdp) );
    tdpfile='';
end


% Load optional arguments
constants = cascadeConstants();
options = constants.defaultMakeplotsOptions;

if numel(args)==2,
    assert( isstruct(args{2}) );
    options = catstruct( options, args{2} );

elseif numel(args)>2,
    args = args(2:end);
    assert( iscell(args) & mod(numel(args),2)==0, ...
            'Incorrect format for optional arguments list' );
    vopt = struct(args{:});
    options = catstruct( options, vopt );
end



f_axis=tdp(2:end,1);
i_axis=tdp(1,2:end);

% Setup contour levels
top = options.tdp_max;
con=0:(top/size(options.cmap,1)):top;
tdp(end,end) = 10;  % hack to make colorscale fixed

% draw contour plot
[~,hand]=contourf(cax,i_axis,f_axis,tdp(2:end,2:end),con);

% Extra formatting
if tdpfile~=0,
    tdpfile = strrep(tdpfile,'_',' ');
    title(tdpfile);
end
set(hand,'LineColor','none');
colormap(cax,options.cmap);
% xlabel(cax,'Initial FRET');
% set(cax,'XTick', i_axis(1):0.1:i_axis(end) )
% set(cax,'YTick', f_axis(1):0.1:f_axis(end) )
set(cax, 'PlotBoxAspectRatio', [1 1 1]);

% if showLHA
%     ylabel(cax,'Final FRET');
% end

zoom on;
