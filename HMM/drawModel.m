function [hFigure,hText,hRate] = drawModel( model, varargin )
%%
%
%



%---- Process input arguments
if nargin<1,
    [fname,pathname] = uigetfile( '*.qmf','Load a QuB model file...');
    if fname==0, return; end
    model = [pathname fname];
end

if ischar(model)
    tree = qub_loadTree(model);
else
    tree = model;
end


%---- Define parameters

% Define colors (same as in QuB).
colors = [0 0 0; 1 0 0; 0 0 1; 0 1 0; 1 1 0; 0.7 0 0.7];
textColors = [1 1 1; 0 0 0; 1 1 1; repmat([0 0 0],[3,1])];



bw = 5; %width of boxes.




%%
% TODO: add context menus to state boxes to 


hFigure = figure;
cla; hold on;



% Define styles for state squares:
textStyle = {'HorizontalAlignment','center','FontWeight','Bold', ...
             'Margin',bw, 'ButtonDownFcn',@showStateDetails};
         
% Draw a box for each state.
states = tree.States.State;
nStates = numel(states);

hText = zeros(nStates,1);
hRect = zeros(nStates,1);
x = zeros(nStates,1);
y = zeros(nStates,1);

for i=1:numel(states),
    x(i) = states(i).x.data;
    y(i) = 90-states(i).y.data;
    class = states(i).Class.data+1;
    
%     hRect(i) = rectangle( 'Position',[x(i) y(i) 0.5 0.5], 'Curvature',[0.25 0.25], ...
%                'FaceColor',colors(class,:));
    
    hText(i) = text( x(i)+bw/2,y(i)+bw/3, num2str(i), ...
          'Color',textColors(class,:), 'BackgroundColor',colors(class,:), ...
          'UserData',i, textStyle{:} );
        
    % Create a context menu for editing state properties
end

axis([0 100 0 100]);


% Draw lines and rate numbers for each 
rates = tree.Rates.Rate;
hRate = zeros( numel(rates), 2 );
% hArrow = [];

for i=1:numel(rates),
    states = rates(i).States.data+1;
    k0  = rates(i).k0.data;
    
    x_d = x(states(1))/2 + ( x(states(2))-x(states(1)) )/4;
    y_d = y(states(1))/2 + ( y(states(2))-y(states(1)) )/4;
    x_arr = x(states)/2+x_d;
    y_arr = y(states)/2+y_d;
    
    arrowOptions = {'Length',8, 'TipAngle',30};
    
    % Determine the angle of the arrow for placement of rate numbers.
    angle = atan( diff(y(states))/diff(x(states)) )*180/pi; %degrees
    angle = abs(mod(angle,180));
    
    % Show rates.
    if angle<45 || angle>135, %horizontal-ish
        arrow( [x_arr(1) y_arr(1)-0.75]+bw/2, [x_arr(2) y_arr(2)-0.75]+bw/2, 'ends','start', arrowOptions{:} );
        arrow( [x_arr(1) y_arr(1)+0.75]+bw/2, [x_arr(2) y_arr(2)+0.75]+bw/2, 'ends','stop',  arrowOptions{:} );
        
        hRate(i,1) = text( mean(x(states))+2,mean(y(states))+4, num2str(k0(1)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom');
        hRate(i,2) = text( mean(x(states))+2,mean(y(states))-1, num2str(k0(2)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','top');
      
    else %vertical-ish
        arrow( [x_arr(1)-0.75 y_arr(1)]+bw/2, [x_arr(2)-0.75 y_arr(2)]+bw/2, 'ends','start', arrowOptions{:} );
        arrow( [x_arr(1)+0.75 y_arr(1)]+bw/2, [x_arr(2)+0.75 y_arr(2)]+bw/2, 'ends','stop',  arrowOptions{:} );
        
        hRate(i,1) = text( mean(x(states))+5,mean(y(states)), num2str(k0(1)), ...
              'HorizontalAlignment','left', 'VerticalAlignment','bottom');
        hRate(i,2) = text( mean(x(states)),mean(y(states)), num2str(k0(2)), ...
              'HorizontalAlignment','right', 'VerticalAlignment','bottom');
    end
    
    % Setup a callback for editing the rates.
    set( hRate(i,:), 'ButtonDownFcn',@editRate );
    set( hRate(i,1), 'UserData',[i 1] );
    set( hRate(i,2), 'UserData',[i 2] );

end


% finish up figure styling
box on;
set( gca, 'xtick',[] );
set( gca, 'ytick',[] );


%%




function showStateDetails( hObject,eventdata,handles )

% Get identifier for rate = index into rates structure in model.
stateID = get(hObject,'UserData');

p0    = tree.States.State(stateID).Pr.data;  %Initial probability
class = tree.States.State(stateID).Class.data+1;
mu    = tree.Amps.data(class); %mean FRET value
sigma = tree.Stds.data(class); %standard deviation of about mean

% Request new rate value from user and update the GUI.
f = inputdlg( {'Initial probability (0..1)','Class number (1..9):',...
               'Mean FRET value:','Stdev of FRET:'}, ...
              'Enter updated state parameters values:',1, ...
              {num2str(p0),num2str(class),num2str(mu),num2str(sigma)});

if isempty(f)
    return;
end

% Update model data
tree.States.State(stateID).Pr.data = str2double(f{1});
tree.States.State(stateID).Class.data = str2double(f{2})-1;
tree.Amps.data(class) = str2double(f{3}); % WARNING these were set for old class!
tree.Stds.data(class) = str2double(f{4});

% If the model has changed, update square colors.
newClass = str2double(f{2});

if class~=newClass,
    set( hText(stateID), 'Color',           textColors(newClass,:), ...
                         'BackgroundColor', colors(newClass,:) );
end

end



function editRate( hObject,eventdata,handles )

    
% Get identifier for rate = index into rates structure in model.
userData = get(hObject,'UserData');
connectionID = userData(1);
direction = userData(2); %index into rate data.

k0 = tree.Rates.Rate(connectionID).k0.data(direction);

% Request new rate value from user and update the GUI.
f = inputdlg('Enter updated rate value (per second):','Specify rate',...
             1,{num2str(k0)});
if ~isempty(f)
    k0 = str2double(f);
    tree.Rates.Rate(connectionID).k0.data(direction) = k0;
    
    % Update GUI
    set( hRate(connectionID,direction), 'String',num2str(k0) );
end


end




end %FUNCTION

