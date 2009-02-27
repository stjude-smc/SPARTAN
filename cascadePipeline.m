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

% Last Modified by GUIDE v2.5 26-Oct-2008 14:57:38

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


% --- Executes on button press in btnGettraces.
function btnGettraces_Callback(hObject, eventdata, handles)
% hObject    handle to btnGettraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gettraces;


% --- Executes on button press in btnAutotrace.
function btnAutotrace_Callback(hObject, eventdata, handles)
% hObject    handle to btnAutotrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
autotrace;


% --- Executes on button press in btnSorttraces.
function btnSorttraces_Callback(hObject, eventdata, handles)
% hObject    handle to btnSorttraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sorttraces;


% --- Executes on button press in btnMakeplots.
function btnMakeplots_Callback(hObject, eventdata, handles)
% hObject    handle to btnMakeplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
makeplots;


% --- Executes on button press in btnOverlay.
function btnOverlay_Callback(hObject, eventdata, handles)
% hObject    handle to btnOverlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frethistComparison;


