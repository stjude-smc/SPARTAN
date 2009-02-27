function tplot(tdp, constants, varargin)
% TPLOT  Draws a transition density contour plot
%
%   TPLOT( TDP, ... )
%   Draws the transition density contour plot in TDP.  If no arguments
%   supplied, the function will prompt the user for the file.
%   

%Plots transition density plot created by tdplot from a file, or from the
%variable tdp entered at command line.

% If not given, load constants
if nargin<2,
    constants = cascadeConstants();
end


% Load TD Plot data
tdpfile=0;
if nargin==0
    [tdpfile tdppath]=uigetfile('_tdp.txt','Choose tdplot:');
    if tdpfile==0
        disp('No File Selected.')
        return
    else
        file=strcat(tdppath,tdpfile);
        tdp=load(file);
    end
end

f_axis=tdp(2:end,1);
i_axis=tdp(1,2:end);

% Load stanford colormap
CMAP_PATH='frethist_colormap.txt';
cmap=dlmread(CMAP_PATH,' ')/255;

% Setup contour levels
% top=max(max(tdp(2:end,2:end)));
% disp(top);
top = constants.tdp_max;
% top = 0.5;
% nl = size(cmap,1)-1;
con=0:(top/13):top;

tdp(end,end) = 1;  % hack to make colorscale fixed

% draw contour plot
% figure;
[C,hand]=contourf(i_axis,f_axis,tdp(2:end,2:end),con);

% Extra formatting
if tdpfile~=0
    title(tdpfile);
end
set(hand,'LineColor','none');
colormap(cmap);
% xlabel('Initial FRET');
% set(gca,'XTick', i_axis(1):0.1:i_axis(end) )
% set(gca,'YTick', f_axis(1):0.1:f_axis(end) )
set(gca, 'PlotBoxAspectRatio', [1 1 1]);

% if showLHA
%     ylabel('Final FRET');
% end

zoom on;
