function varargout = spartan(varargin)
% SPARTAN M-file for spartan.fig
%      SPARTAN, by itself, creates a new SPARTAN or raises the existing
%      singleton*.
%
%      H = SPARTAN returns the handle to a new SPARTAN or the handle to
%      the existing singleton*.
%
%      SPARTAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPARTAN.M with the given input arguments.
%
%      SPARTAN('Property','Value',...) creates a new SPARTAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spartan_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spartan_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% Last Modified by GUIDE v2.5 14-Aug-2015 16:35:04


% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spartan_OpeningFcn, ...
                   'gui_OutputFcn',  @spartan_OutputFcn, ...
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


function listPrograms() %#ok<DEFNU>
% This function is never called. Its purpose is to list all functions used by
% the program, which are actually in the callbacks in the .fig file.
% This is to make compiling to program easier.
gettraces; autotrace; rtdgui;
sorttraces; makeplots; frethistComparison; avgFretTime; tsqComparison;
batchKinetics; simulate; lifetime_exp; percentTime;
crosstalkcorrect; gammacorrect; scaleacceptor; haranfilter;
forQuB; fluorToFret; forHammy; forvbFRET; hammyToDWT; vbFRET_dwt; forOrigin;
combineDatasets; resizeTraces;




% --- Executes just before spartan is made visible.
function spartan_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spartan (see VARARGIN)

constants = cascadeConstants;
set( handles.figure1, 'Name', constants.software );
updateSpartan; %check for updates

% Set working directory in GUI
set(handles.txtCWD, 'String',pwd);

% Choose default command line output for spartan
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spartan wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spartan_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




function txtCWD_Callback(hObject, ~, handles) %#ok<DEFNU>
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
function btnBrowse_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to btnBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = uigetdir('Select working directory');
if d==0, return; end

cd(d);
set(handles.txtCWD, 'String',pwd);

guidata(hObject, handles);
