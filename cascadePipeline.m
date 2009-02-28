function varargout = cascadePipeline(varargin)
% CASCADEPIPELINE M-file for cascadePipeline.fig
%      CASCADEPIPELINE, by itself, creates a new CASCADEPIPELINE or raises the existing
%      singleton*.
%
%      H = CASCADEPIPELINE returns the handle to a new CASCADEPIPELINE or the handle to
%      the existing singleton*.
%
%      CASCADEPIPELINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CASCADEPIPELINE.M with the given input arguments.
%
%      CASCADEPIPELINE('Property','Value',...) creates a new CASCADEPIPELINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cascadePipeline_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cascadePipeline_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cascadePipeline

% Last Modified by GUIDE v2.5 28-Feb-2009 14:53:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cascadePipeline_OpeningFcn, ...
                   'gui_OutputFcn',  @cascadePipeline_OutputFcn, ...
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


% --- Executes just before cascadePipeline is made visible.
function cascadePipeline_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cascadePipeline (see VARARGIN)

% Set working directory in GUI
set(handles.txtCWD, 'String',pwd);

% Choose default command line output for cascadePipeline
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cascadePipeline wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cascadePipeline_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




function txtCWD_Callback(hObject, eventdata, handles)
% hObject    handle to txtCWD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtCWD as text
%        str2double(get(hObject,'String')) returns contents of txtCWD as a double

d = get(hObject,'String');
if ~exist(d,'dir')
    warning('Directory does not exist!');
    set(handles.txtCWD, 'String',pwd);
else
    cd(d);
end



% --- Executes on button press in btnBrowse.
function btnBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to btnBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = uigetdir('Select working directory');
cd(d);
set(handles.txtCWD, 'String',pwd);

guidata(hObject, handles);


% --- Executes on button press in btnKinetics.
function btnKinetics_Callback(hObject, eventdata, handles)
% hObject    handle to btnKinetics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnSimulate.
function btnSimulate_Callback(hObject, eventdata, handles)
% hObject    handle to btnSimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnOverlay2.
function btnOverlay2_Callback(hObject, eventdata, handles)
% hObject    handle to btnOverlay2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


