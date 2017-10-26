function output = settingdlg(input, varargin) %fields, Prompt, types, Title, fun, params)
%settingdlg  Input dialog to change arbitrary struct fields
%
%   OUT = settingdlg(IN) prompts the user for values in the struct IN,
%   saving them in OUT if the user hits "OK" and returns IN unmodified otherwise.
%   User input is converted into the type given in the 
%
%   OUT = settingdlg(IN,FIELDS) specifies which fields in OPT to show.
%   Other fields will be passed to OUT unmodified.
%
%   OUT = settingdlg(IN,FIELDS,PROMPTS) specifies a title for each
%   input field.
% 
%   OUT = settingdlg(IN,FIELDS,PROMPTS,TYPES) specifies valid input
%   values for each field. Each element of the TYPES cell array can be
%   empty (use input type), a cell array of strings (), or a handle to a
%   function that returns true if the input is valid and false otherwise.
%
%   Example:
%     in     = struct('num',11, 'log',false, 'str','a');
%     fields = fieldnames(in);  %all fields
%     prompt = {'Number','Logical','Choice'};
%     types  = { @isfinite, [], {'a','b','c'} };
%     out    = settingdlg(in, fields, prompt, types, 'Settings');
%
%   settingdlg(...,FUN,FIN) after the dialog closes, calls the function
%   handle as FUN(OPT,FIN{:}).
%
%   See also: inputdlg.

%   Based on MATLAB's built-in inputdlg function.
%   Copyright 2015-2016 Cornell University All Rights Reserved.



%% Process input arguments
narginchk(1,8);
nargoutchk(0,1);

if ~isstruct(input) || numel(input)>1,
    error('First argument must be a scalar struct');
end

% Assign default values
fields = fieldnames(input);
Prompt = strcat(fields,':');
types  = cell(numel(fields),1);
Title  = 'Settings';

% Parse out field and function arguments from beginning and end of varargin.
idxFun = find(  cellfun( @(x)isa(x,'function_handle'), varargin )  );
assert( numel(idxFun)<=1, 'Multiple function inputs not allowed' );
if ~isempty(idxFun)
    [fun,fin] = varargin{idxFun:end};
    varargin = varargin(1:idxFun-1);
end

if numel(varargin)>=1,
    fields = varargin{1};
end
if numel(varargin)>=2,
    Prompt = varargin{2};
end
if numel(varargin)>=3,
    types  = varargin{3};
end
if numel(varargin)>=4,
    Title  = varargin{4};
end


NumQuest = numel(Prompt);

if numel(types) < NumQuest,
    [types{numel(types)+1:NumQuest}] = deal([]);  %fill trailing missing values
end



%% Create Figure and calculate control properties

% Create figure
FigWidth = 175;
FigColor = get(0,'DefaultUicontrolBackgroundColor');

InputFig = dialog(                 ...
  'Visible'          ,'off'      , ...  %wait until fully assembled
  'KeyPressFcn'      ,@doFigureKeyPress, ...
  'Name'             ,Title      , ...
  'Pointer'          ,'arrow'    , ...
  'Units'            ,'pixels'   , ...
  'UserData'         ,'Cancel'   , ...  %cancel by default if window closed
  'Tag'              ,Title      , ...
  'HandleVisibility' ,'callback' , ...
  'Color'            ,FigColor   , ...
  'NextPlot'         ,'add'      , ...
  'WindowStyle'      ,'modal'    , ...
  'Resize'           ,'off'       ...
  );

% Set default control properties
DefOffset    = 5;
DefBtnWidth  = 53;
DefBtnHeight = 23;

TextInfo.Units              = 'pixels'   ;
TextInfo.FontSize           = get(0,'FactoryUicontrolFontSize');
TextInfo.FontWeight         = get(InputFig,'DefaultTextFontWeight');
TextInfo.HorizontalAlignment= 'left'     ;
TextInfo.HandleVisibility   = 'callback' ;

StInfo=TextInfo;
StInfo.Style              = 'text'  ;
StInfo.BackgroundColor    = FigColor;

EdInfo=StInfo;
EdInfo.BackgroundColor = 'white';

BtnInfo=StInfo;
BtnInfo.Style               = 'pushbutton';
BtnInfo.HorizontalAlignment = 'center';
BtnInfo.KeyPressFcn         = @doControlKeyPress;
BtnInfo.Callback            = @doCallback;

% Add VerticalAlignment here as it is not applicable to the above.
TextInfo.VerticalAlignment  = 'bottom';
TextInfo.Color              = get(0,'FactoryUicontrolForegroundColor');



%% Create a dummy control to get default control positioning/size.

% Button size
btnMargin=1.4;
ExtControl=uicontrol(InputFig   ,BtnInfo     , ...
  'String'   ,getString(message('MATLAB:uistring:popupdialogs:Cancel'))        , ...
  'Visible'  ,'off'         ...
  );

BtnExtent = get(ExtControl,'Extent');
BtnWidth  = max(DefBtnWidth,BtnExtent(3)+8);
BtnHeight = max(DefBtnHeight,BtnExtent(4)*btnMargin);
delete(ExtControl);


% Textbox size, taking into account multi-line (word-wrapped) prompts
TxtWidth=FigWidth-2*DefOffset;
ExtControl = uicontrol(InputFig, StInfo, 'String','', ...
  'Position' ,[ DefOffset DefOffset 0.96*TxtWidth BtnHeight ] , ...
  'Visible'  ,'off'        ...
  );

WrapQuest=cell(NumQuest,1);
QuestPos=zeros(NumQuest,4);

for ExtLp=1:NumQuest
    [WrapQuest{ExtLp},QuestPos(ExtLp,:)]= ...
      textwrap(ExtControl,Prompt(ExtLp),80);
end % for ExtLp
delete(ExtControl);

QuestWidth =QuestPos(:,3);
QuestHeight=QuestPos(:,4);

if ismac % Change Edit box height to avoid clipping on mac.
    editBoxHeightScalingFactor = 1.4;
else 
    editBoxHeightScalingFactor = 1;
end
TxtHeight = QuestHeight(1)/size(WrapQuest{1,1},1) * editBoxHeightScalingFactor;
EditHeight = TxtHeight+4;

FigHeight = (NumQuest+2)*DefOffset + BtnHeight+(NumQuest*EditHeight) + sum(QuestHeight);

QuestYOffset = zeros(NumQuest,1);
EditYOffset  = zeros(NumQuest,1);
QuestYOffset(1) = FigHeight-DefOffset-QuestHeight(1);
EditYOffset(1)  = QuestYOffset(1)-EditHeight;

for YOffLp=2:NumQuest,
    QuestYOffset(YOffLp) = EditYOffset(YOffLp-1) - QuestHeight(YOffLp) - DefOffset;
    EditYOffset(YOffLp)  = QuestYOffset(YOffLp)  - EditHeight;
end % for YOffLp



%% Create input controls
AxesHandle  = axes('Parent',InputFig, 'Position',[0 0 1 1], 'Visible','off');
EditHandle  = zeros(NumQuest,1);
QuestHandle = zeros(NumQuest,1);

for lp=1:NumQuest,
    DefAns = input.(fields{lp});
    if isnumeric(DefAns),
        DefAns = num2str(DefAns);
    end
  
    if iscell(types{lp}) && ~isempty(types{lp}),
        val = find(strcmpi(DefAns,types{lp}));
        if isempty(val),
            warning('Invalid default value for field %s',fields{lp});
            val = 1;
        end

        EditHandle(lp) = uicontrol(InputFig, EdInfo, 'Style','popup', ...
            'String'  ,types{lp}, ...
            'Value'   ,val );
        
    elseif islogical(DefAns)
        EditHandle(lp) = uicontrol(InputFig, EdInfo, 'Style','checkbox', ...
            'BackgroundColor', FigColor, ...
            'Value'      , DefAns );
    else
        EditHandle(lp) = uicontrol(InputFig, EdInfo, 'Style','edit', ...
            'String'     ,DefAns );
    end
    set( EditHandle(lp), 'Position',[DefOffset EditYOffset(lp) TxtWidth EditHeight] );

    QuestHandle(lp) = text('Parent',AxesHandle, TextInfo, ...
        'Position'   ,[DefOffset QuestYOffset(lp)], ...
        'String'     ,WrapQuest{lp} );

    MinWidth = max(QuestWidth(:));
    
    % Get the extent of the text object. See g1008152
    questExtent = get(QuestHandle(lp), 'Extent');
    MinWidth = max(MinWidth, questExtent(3));  
    FigWidth = max(FigWidth, MinWidth+2*DefOffset);

end % for lp

% fig width may have changed, update the edit fields
TxtWidth = FigWidth-2*DefOffset;
for lp=1:NumQuest
    set(EditHandle(lp), 'Position', [DefOffset EditYOffset(lp) TxtWidth EditHeight]);
end

FigWidth = max(FigWidth,2*(BtnWidth+DefOffset)+DefOffset);
FigPos = [0 0 FigWidth FigHeight];
set(InputFig,'Position',getnicedialoglocation(FigPos,get(InputFig,'Units')));

uicontrol(InputFig, BtnInfo, ...
  'Position'   ,[ FigWidth-2*BtnWidth-2*DefOffset DefOffset BtnWidth BtnHeight ] , ...
  'String'     ,'OK', ...
  'Tag'        ,'OK'        , ...
  'UserData'   ,'OK'          ...
  );

% setdefaultbutton2(InputFig, OKHandle);  %highlight default button

uicontrol(InputFig, BtnInfo, ...
  'Position'   ,[ FigWidth-BtnWidth-DefOffset DefOffset BtnWidth BtnHeight ]           , ...
  'String'     ,'Cancel', ...
  'Tag'        ,'Cancel'    , ...
  'UserData'   ,'Cancel'      ...
  );



%% Show dialog and wait for user input.
movegui(InputFig);  %make sure dialog is entirely on screen
set(InputFig,'Visible','on');
drawnow;

% Set focus to first input control
if ~isempty(EditHandle)
  uicontrol(EditHandle(1));
end



%% Wait for user input, only accepting 'OK' if all entries are valid.
% Cancel returns the input. Should it be empty?
output = input;
value = cell(size(types));
ok = false;

while ishghandle(InputFig)
    % Wait for user to close the window or press a button.
    uiwait(InputFig);
    
    % Close dialog with no output if user closes dialog or clicks Cancel.
    if ~ishghandle(InputFig), break; end
    if strcmp(get(InputFig,'UserData'),'Cancel'), break; end

    % User clicked OK. Verify input is correct. If not, try again.
    % FIXME: handle fields with multiple allowed types (number, string).
    valid = true(size(types));
    
    for i=1:NumQuest,
        DefAns = input.(fields{i});
        
        try
            % Get string value from dropdown combo boxes.
            switch lower(get(EditHandle(i),'Style'))
                case 'popupmenu'
                    value{i} = types{i}{ get(EditHandle(i),'Value') };
                case 'checkbox'
                    value{i} = logical( get(EditHandle(i),'Value') );
                otherwise
                    value{i} = get(EditHandle(i),'String');
            end

            % Convert to numeric inputs to double
            if isnumeric(DefAns) && ischar(value{i})
                value{i} = str2double(value{i});
                valid(i) = ~isnan(value{i});
            end

            % Check validity of the final input
            if isa(types{i},'function_handle'),
                valid(i) = types{i}(value{i});
            end
        catch
            valid(i) = false;
        end
    end %for each dialog field
    

    if all(valid),
        % All inputs are valid. Exit returning inputs.
        for i=1:NumQuest,
            output.(fields{i}) = value{i};
        end
        delete(InputFig);
        ok = true;
        
    else
        % Highlight invalid fields and prompt user to fix them
        errordlg('Highlighted fields are invalid.','Invalid value');
        set( QuestHandle(valid),  'Color','k' );
        set( QuestHandle(~valid), 'Color','r' );
    end
end
drawnow; % Update the view to remove the closed figure (g1031998)


% Execute callback to update target with new settings, if requested.
if ok && ~isempty(idxFun),
    fun(fin{:}, output);
end
if ~ok, output=[]; end


end



%% Callback functions
function doFigureKeyPress(~, evd)
switch(evd.Key)
  case {'return','space'}
    set(gcbf,'UserData','OK');
    uiresume(gcbf);
  case {'escape'}
    delete(gcbf);
end
end

function doControlKeyPress(obj, evd)
switch(evd.Key)
  case {'return'}
    if ~strcmp(get(obj,'UserData'),'Cancel')
      set(gcbf,'UserData','OK');
      uiresume(gcbf);
    else
      delete(gcbf)
    end
  case 'escape'
    delete(gcbf)
end
end

function doCallback(obj, ~)
if ~strcmp(get(obj,'UserData'),'Cancel')
  set(gcbf,'UserData','OK');
  uiresume(gcbf);
else
  delete(gcbf)
end
end



%% Accessory functions

function figure_size = getnicedialoglocation(figure_size, figure_units)
% adjust the specified figure position to fig nicely over GCBF
% or into the upper 3rd of the screen

%  Copyright 1999-2011 The MathWorks, Inc.

parentHandle = gcbf;
convertData.destinationUnits = figure_units;
convertData.reference = 0;
if ~isempty(parentHandle)
    % If there is a parent figure
    convertData.hFig = parentHandle;
    convertData.size = get(parentHandle,'Position');
    convertData.sourceUnits = get(parentHandle,'Units');  
    c = []; 
else
    % If there is no parent figure, use the root's data
    % and create a invisible figure as parent
    convertData.hFig = figure('visible','off');
    convertData.size = get(0,'ScreenSize');
    convertData.sourceUnits = get(0,'Units');
    c = onCleanup(@() close(convertData.hFig));
end

% Get the size of the dialog parent in the dialog units
container_size = hgconvertunits(convertData.hFig, convertData.size ,...
    convertData.sourceUnits, convertData.destinationUnits, convertData.reference);

delete(c);

figure_size(1) = container_size(1)  + 1/2*(container_size(3) - figure_size(3));
figure_size(2) = container_size(2)  + 2/3*(container_size(4) - figure_size(4));

end





