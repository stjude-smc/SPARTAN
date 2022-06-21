function [ax,AxisPos] = subplotTight(nRow, nCol, varargin)
% subplotTight  faster version of subplot with tight plot spacing
%

% https://www.mathworks.com/matlabcentral/answers/16458-making-less-space-between-figures-in-subplot


%% Process input arguments
narginchk(2,Inf);
nargoutchk(0,1);

% Spacing between rows: leave room for axis ticks and labels.
BplusT = 0.1 + 0.05 * nRow;

% Spacing between columns: leave no room (label and ticks hidden).
LplusR = 0.02;

% Spacing from edges (drawable area): [left bottom width height]
Lplus = (sqrt(0.12)/nCol)^2;
defPos = [0.1+Lplus, 0.1, 0.845-Lplus, 0.845];


%% Define axis positions
nPlot = nRow * nCol;
plots = 0:(nPlot - 1);
row   = (nRow - 1) - fix(plots(:) / nCol);
col   = rem(plots(:), nCol);
col_offset  = defPos(3) * LplusR / (nCol - LplusR);
row_offset  = defPos(4) * BplusT / (nRow - BplusT);
totalwidth  = defPos(3) + col_offset;
totalheight = defPos(4) + row_offset;
width       = totalwidth  / nCol - col_offset;
height      = totalheight / nRow - row_offset;

if width * 2 > totalwidth / nCol
   if height * 2 > totalheight / nRow
      AxisPos = [(defPos(1) + col * totalwidth / nCol), ...
            (defPos(2) + row * totalheight / nRow), ...
            width(ones(nPlot, 1), 1), ...
            height(ones(nPlot, 1), 1)];
   else
       AxisPos = [(defPos(1) + col * totalwidth / nCol), ...
            (defPos(2) + row * defPos(4) / nRow), ...
            width(ones(nPlot, 1), 1), ...
            (0.7 * defPos(ones(nPlot, 1), 4) / nRow)];
   end
else
   if height * 2 <= totalheight / nRow
      AxisPos = [(defPos(1) + col * defPos(3) / nCol), ...
            (defPos(2) + row * defPos(4) / nRow), ...
            (0.7 * defPos(ones(nPlot, 1), 3) / nCol), ...
            (0.7 * defPos(ones(nPlot, 1), 4) / nRow)];
   else
      AxisPos = [(defPos(1) + col * defPos(3) / nCol), ...
            (defPos(2) + row * totalheight / nRow), ...
            (0.7 * defPos(ones(nPlot, 1), 3) / nCol), ...
            height(ones(nPlot, 1), 1)];
    end
end


%% Create axes

ax = zeros( nCol, nRow );
for i = 1:nPlot
  ax(i) = axes( 'Position',AxisPos(i, :), 'Box','on', varargin{:} );
end
ax = ax';


end  %function



