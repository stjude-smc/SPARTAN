function varargout = rtdgui(varargin)
% RTDGUI MATLAB code for rtdgui.fig
%      RTDGUI, by itself, creates a new RTDGUI or raises the existing
%      singleton*.
%
%      H = RTDGUI returns the handle to a new RTDGUI or the handle to
%      the existing singleton*.
%
%      RTDGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RTDGUI.M with the given input arguments.
%
%      RTDGUI('Property','Value',...) creates a new RTDGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rtdgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rtdgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rtdgui

% Last Modified by GUIDE v2.5 31-Aug-2015 13:04:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rtdgui_OpeningFcn, ...
                   'gui_OutputFcn',  @rtdgui_OutputFcn, ...
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


% --- Executes just before rtdgui is made visible.
function rtdgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rtdgui (see VARARGIN)

% Choose default command line output for rtdgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
constants = cascadeConstants;

% set default values
% FIXME: these and the defaults in rtdTool should come from some shared
% location so that they can all be changed from one place.
modelfile = fullfile(constants.modelLocation,'tRNA selection','2014_04_18 EColi.qmf');
if ~exist(modelfile,'file'),
    [~,f,e] = fileparts(modelfile);
    modelfile = [f,e];
end
set(handles.editModelPath,'String',modelfile);

set(handles.selProdDwell,'String','120');
updateStateDropdown(hObject,eventdata,handles);

set(handles.selProdState,'Value',4);
set(handles.eventFilter,'Value',1);
set(handles.selProdState,'Enable','on');
set(handles.selProdDwell,'Enable','on');

set(handles.advancedSettings,'String','');
set(handles.useAdvSettings,'Value',0);
set(handles.advancedSettings,'Enable','off');
set(handles.browseSettings,'Enable','off');

% Can't load .m files in compiled code, so this feature is disabled.
if isdeployed
    set(handles.advancedSettings,'Enable','off');
end

% UIWAIT makes rtdgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function updateStateDropdown(hObject,eventdata,handles)
modelPath = get(handles.editModelPath,'String');
if exist(modelPath,'file')==2
    qubModel = qub_loadModel(modelPath);
    for i=1:qubModel.nClass
        stateList{i} = num2str(i);
    end
    set(handles.selProdState,'Value',1);
    set(handles.selProdState,'String',stateList);
else
    errordlg('Invalid model file.');
end

% --- Outputs from this function are returned to the command line.
function varargout = rtdgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editModelPath_Callback(hObject, eventdata, handles)
% hObject    handle to editModelPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editModelPath as text
%        str2double(get(hObject,'String')) returns contents of editModelPath as a double
updateStateDropdown(hObject,eventdata,handles)


% --- Executes on button press in browseModelPath.
function browseModelPath_Callback(hObject, eventdata, handles)
% hObject    handle to browseModelPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldPath = get(handles.editModelPath,'String');
[FileName,PathName,FilterIndex] = uigetfile({'*.qmf','QuB model file'},'Select a model file for idealization.',oldPath);
if FilterIndex
    set(handles.editModelPath,'String',fullfile(PathName,FileName));
    updateStateDropdown(hObject,eventdata,handles);
end


% --- Executes on button press in eventFilter.
function eventFilter_Callback(hObject, eventdata, handles)
% hObject    handle to eventFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eventFilter

filterOn = get(handles.eventFilter,'Value');
if filterOn
    set(handles.selProdState,'Enable','on');
    set(handles.selProdDwell,'Enable','on');
else
    set(handles.selProdState,'Enable','off');
    set(handles.selProdDwell,'Enable','off');
end   


% --- Executes on selection change in selProdState.
function selProdState_Callback(hObject, eventdata, handles)
% hObject    handle to selProdState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selProdState contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selProdState



function selProdDwell_Callback(hObject, eventdata, handles)
% hObject    handle to selProdDwell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of selProdDwell as text
%        str2double(get(hObject,'String')) returns contents of selProdDwell as a double



% --- Executes on button press in okBtn.
function okBtn_Callback(hObject, eventdata, handles)
% hObject    handle to okBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
opt.kinModel = get(handles.editModelPath,'String');
opt.prodDwell = str2num(get(handles.selProdDwell,'String'));
opt.prodState = get(handles.selProdState,'Value');
opt.skipStateFilter = ~get(handles.eventFilter,'Value');
if get(handles.useAdvSettings,'Value') && ~isdeployed
    settingsFile = get(handles.advancedSettings,'String');
    run(settingsFile);
    if exist('customOpt','var')
        opt = mergestruct(opt,customOpt);
    else
        errordlg('Advanced settings file does not contain valid customOpt struct. Proceeding with default settings.');
    end
end
rtdTool(opt);




% --- Executes on button press in cancelBtn.
function cancelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to cancelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);


% --- Executes on button press in useAdvSettings.
function useAdvSettings_Callback(hObject, eventdata, handles)
% hObject    handle to useAdvSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useAdvSettings
if get(handles.useAdvSettings,'Value') && ~isdeployed;
    set(handles.advancedSettings,'Enable','on');
    set(handles.browseSettings,'Enable','on');
else
    set(handles.advancedSettings,'Enable','off');
    set(handles.browseSettings,'Enable','off');
end



function advancedSettings_Callback(hObject, eventdata, handles)
% hObject    handle to advancedSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of advancedSettings as text
%        str2double(get(hObject,'String')) returns contents of advancedSettings as a double



% --- Executes on button press in browseSettings.
function browseSettings_Callback(hObject, eventdata, handles)
% hObject    handle to browseSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldPath = get(handles.advancedSettings,'String');
[FileName,PathName,FilterIndex] = uigetfile({'*.m','Advanced settings script'},'Select the MATLAB script to use for advanced settings.',oldPath);
if FilterIndex
    set(handles.advancedSettings,'String',fullfile(PathName,FileName));
end


% --- Executes on button press in remakePlots.
function remakePlots_Callback(hObject, eventdata, handles)
% hObject    handle to remakePlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
opt = struct;
if get(handles.useAdvSettings,'Value')
    settingsFile = get(handles.advancedSettings,'String');
    run(settingsFile);
    if exist('customOpt','var')
        opt = mergestruct(opt,customOpt);
    else
        errordlg('Advanced settings file does not contain valid customOpt struct. Proceeding with default settings.');
    end
end
rtdPlots(opt);
