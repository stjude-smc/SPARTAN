function varargout = qubstatedlg(varargin)
% QUBSTATEDLG  Dialog to change parameters of a QubModel object
%
%   MODEL = qubstatedlg(MODEL, STATE) opens a modal dialog to modify model 
%   parameters of the state indicated by the integer STATE of the QubModel 
%   object MODEL.
% 
%   See also: QubModelViewer, QubModel, batchKinetis.

%   Copyright 2018 Cornell University All Rights Reserved.

% Last Modified by GUIDE v2.5 19-Sep-2018 13:52:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @qubstatedlg_OpeningFcn, ...
                   'gui_OutputFcn',  @qubstatedlg_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before qubstatedlg is made visible.
function qubstatedlg_OpeningFcn(hObject, ~, handles, model, state)
% 

handles.output = 'Cancel';    %default dialog result if window closed

handles.model = model;
guidata(hObject, handles);

% Set dialog position on screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    % Center dialog on screen
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);
    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    % Center dialog over parent window
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

set(handles.figure1,'WindowStyle','modal')

% % Preview state number & color
% set( handles.axes1, ...
%     'Visible','off', 'YDir','reverse', ...
%     'XLim', get(Img,'XData'),  'YLim', get(Img,'YData')  );

% Set dialog control values here from inputModel.
classList = cellfun(@num2str, num2cell(1:model.nClasses+1), 'Uniform',false);
set( handles.cboClass, 'String',classList, 'Value',model.class(state) );

set( handles.edStateNo,  'String', num2str(state) );
set( handles.edStartProb,'String', num2str(model.p0(state)) );

cboClass_Callback([], [], handles);

% Wait for user input
uiwait(handles.figure1);

% END FUNCTION qubstatedlg_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = qubstatedlg_OutputFcn(~, ~, handles)
% Return updated QubModel to user, or empty if "Cancel" is pressed.

if ~strcmpi(handles.output,'Ok')
    [varargout{1:nargout}] = deal([]);
    delete(handles.figure1);
    return;
end

varargout{1} = handles.model;
delete(handles.figure1);

% END FUNCTION qubstatedlg_OutputFcn



% --- Executes when either 'ok' or 'cancel buttons are pressed.
function btnOk_Callback(hObject, ~, handles)
% 

if strcmpi( get(hObject,'String'), 'Ok' )

    % Get GUI control values
    state = str2double( get(handles.edStateNo,'String') );
    class = get(handles.cboClass,'Value');
    handles.model.class(state) = class;

    p0    = str2double( get(handles.edStartProb,'String') );
    mu    = str2double( get(handles.edMeanFret, 'String') );
    sigma = str2double( get(handles.EdStdFret,  'String') );

    % Check if input parameter values are valid
    valid(1) = ~isnan(p0) && p0>=0 && p0<=1;
    valid(2) = ~isnan(mu) && ~isinf(mu);
    valid(3) = ~isnan(sigma) && ~isinf(sigma) && sigma>=0;

    controls = [handles.edStartProb handles.edMeanFret handles.EdStdFret];
    for i=1:numel(valid)
        set( controls(i), 'ForegroundColor',~valid(i)*[1 0 0] );
    end
    if any(~valid)
        errordlg('Highlighted fields are not valid');
        return;
    end

    % Save results
    handles.model.p0(state)    = p0;
    handles.model.mu(class)    = mu;
    handles.model.sigma(class) = sigma;

    handles.model.fixMu(class)    = get( handles.chkFixMu,   'Value' );
    handles.model.fixSigma(class) = get( handles.chkFixSigma,'Value' );
end

% Close the dialog and return current parameters
handles.output = get(hObject,'String');
guidata(hObject, handles);
uiresume(handles.figure1);

% END FUNCTION btnOk_Callback



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, ~) %#ok<DEFNU>
% 
if isequal(get(hObject, 'waitstatus'), 'waiting')
    uiresume(hObject);
else
    delete(hObject);
end



% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, ~, handles) %#ok<DEFNU>
% 
switch get(hObject,'CurrentKey')
    case 'escape'
        guidata(hObject, handles);
        uiresume(handles.figure1);
    case 'return'
        btnOk_Callback(hObject,[],handles);
end    



% --- Executes on selection change in cboClass.
function cboClass_Callback(~, ~, handles)
% Load model parameters from the new class 
class = get(handles.cboClass, 'Value');

if class > handles.model.nClasses
    handles.model.addClass( handles.model.nClasses+1 );
end
set( handles.edMeanFret, 'String', num2str(handles.model.mu(class))    );
set( handles.EdStdFret,  'String', num2str(handles.model.sigma(class)) );

set( handles.chkFixMu,   'Value', handles.model.fixMu(class)   );
set( handles.chkFixSigma,'Value', handles.model.fixSigma(class) );

% END FUNCTION cboClass_Callback
