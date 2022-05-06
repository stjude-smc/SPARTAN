function varargout = fieldArrangementDialog(varargin)
% FIELDARRANGEMENTDIALOG MATLAB code for fieldArrangementDialog.fig
%      FIELDARRANGEMENTDIALOG, by itself, creates a new FIELDARRANGEMENTDIALOG or raises the existing
%      singleton*.
%
%      H = FIELDARRANGEMENTDIALOG returns the handle to a new FIELDARRANGEMENTDIALOG or the handle to
%      the existing singleton*.
%
%      FIELDARRANGEMENTDIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIELDARRANGEMENTDIALOG.M with the given input arguments.
%
%      FIELDARRANGEMENTDIALOG('Property','Value',...) creates a new FIELDARRANGEMENTDIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fieldArrangementDialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fieldArrangementDialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fieldArrangementDialog

% Last Modified by GUIDE v2.5 05-May-2022 17:43:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fieldArrangementDialog_OpeningFcn, ...
                   'gui_OutputFcn',  @fieldArrangementDialog_OutputFcn, ...
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


% --- Executes just before fieldArrangementDialog is made visible.
function fieldArrangementDialog_OpeningFcn(hObject, ~, handles, varargin)
% 

% Process input arguments
narginchk(5,5);
geo = varargin{1};  %field arrangement matrix with indexes into validFields below
validFields = varargin{2};  %all allowed channel names

% Set default command line output
handles.output = [];

% Handles to drop-down controls to show for each possible field organization.
% 1x1, 1x2, 2x1, 2x2, Concatinated Frames (2,3,4)
handles.toShow = { [handles.cboField11], ...
                   [handles.cboField11 handles.cboField12], ...
                   [handles.cboField11 handles.cboField21], ...
                   [handles.cboField11 handles.cboField12 handles.cboField21 handles.cboField22], ...
                   [handles.cboField11 handles.cboField21], ...
                   [handles.cboField11 handles.cboField21 handles.cboField31], ...
                   [handles.cboField11 handles.cboField21 handles.cboField31 handles.cboField41], ...
                 };

% Fill all cbFields with valid choices
handles.hAll = [handles.cboField11 handles.cboField12 handles.cboField21 handles.cboField22 ...
                handles.cboField31 handles.cboField41];
set( handles.hAll, 'String',[{''} validFields] );
set( handles.hAll, 'Visible','off', 'Value',1 );

% Set field name dropdowns: channels are saved as concatinated frames
if size(geo,3)>1
    assert( size(geo,3)<=4, 'Up to 4 channels supported as concatinated frames.' );
    assert( size(geo,1)==1&size(geo,2)==1, 'Fields may be tiled or concatinated but not both.' );
    
    hCbo = [handles.cboField11 handles.cboField21 handles.cboField31 handles.cboField41];
    for i=1:numel(geo)
        set( hCbo(i), 'Value',geo(i)+1, 'Visible','on' );
    end
    
% Set field name dropdowns: channels are tiled within each frame
else
    hCbo = [handles.cboField11 handles.cboField12 handles.cboField21 handles.cboField22];
    for i=1:numel(hCbo)
        id = get( hCbo(i), 'UserData' );  %index of target element in geo
        if all(id<=size(geo))
            set( hCbo(i), 'Value',geo(id(1),id(2))+1, 'Visible','on' );
        end
    end
end

% Set initial field organization choice
if numel(geo)==1
    sel=1;
elseif size(geo,3)>1
    sel = size(geo,3)+3;
elseif size(geo,1)==1 && size(geo,2)==2
    sel=2;
elseif size(geo,1)==2 && size(geo,2)==1
    sel=3;
elseif size(geo,1)==2 && size(geo,2)==2
    sel=4;
end
set( handles.cboFieldOrganization, 'Value',sel );

% UIWAIT makes fieldArrangementDialog wait for user response (see UIRESUME)
guidata(hObject, handles);
uiwait(handles.figure1);

%end function


% --- Outputs from this function are returned to the command line.
function varargout = fieldArrangementDialog_OutputFcn(~, ~, handles) 
% 

if ~strcmpi(handles.output,'Ok')
    [varargout{1:nargout}] = deal([]);
    delete(handles.figure1);
    return;
end

% Assemble output
geo = [];

for i=1:numel(handles.hAll)
    hCbo = handles.hAll(i);
    id = get(hCbo,'UserData');  %row, col coordinate stored here
    value = get(hCbo,'Value');
    if strcmpi(get(hCbo,'Visible'),'on')
        geo( id(1), id(2) ) = value-1;
    end
end

if get(handles.cboFieldOrganization,'Value')>=5
    geo = reshape(geo, [1 1 numel(geo)]);
end
varargout{1} = geo;

sel = get( handles.hAll, 'Value' );
sel = [sel{:}];
sel = sel(sel>1);
validNames = get(handles.cboField11,'String');
varargout{2} = validNames( sort(sel) );  %automatically in wavelength order.

% [varargout{1:nargout}] = handles.output{1:nargout};
delete(handles.figure1);

%end function


% --- Executes on button press in btnOk.
function btnOk_Callback(hObject, ~, handles)
% Close the dialog and return current parameters

% Verify current state
sel = get( handles.hAll, 'Value' );
sel = [sel{:}];
sel = sel(sel>1);
if numel(sel)>numel(unique(sel))
    errordlg('The same field cannot be selected more than once!');
    return;
end

% Otherwise, return to user.
handles.output = get(hObject,'String');
guidata(hObject, handles);
uiresume(handles.figure1);

% end function


% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, ~, handles)
% Close the dialog and return current parameters
handles.output = get(hObject,'String');
guidata(hObject, handles);
uiresume(handles.figure1);
%end function



% --- Executes on selection change in cboFieldOrganization.
function cboFieldOrganization_Callback(hObject, ~, handles) %#ok<*DEFNU>
% Reset, showing only drop-down boxes associated with selected field
% organization.
sel = get(hObject,'Value');
set( handles.hAll, 'Visible','off', 'Value',1 );
set( handles.toShow{sel}, 'Visible','on' );

%end function cboFieldOrganization_Callback



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, ~)
% 
if isequal(get(hObject, 'waitstatus'), 'waiting')
    uiresume(hObject);
else
    delete(hObject);
end



% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, ~, handles)
% 
switch get(hObject,'CurrentKey')
    case 'escape'
        guidata(hObject, handles);
        uiresume(handles.figure1);
    case 'return'
        btnOk_Callback(hObject,[],handles);
end
