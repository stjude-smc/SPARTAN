function varargout = realtimeAnalysis_settings(varargin)
% REALTIMEANALYSIS_SETTINGS  Trace processing and filtering
%
%   Loads fluorescence traces files (produced using gettraces), makes
%   corrections for crosstalk, background intensity, and sets FRET=0
%   where the donor dye is dark.  Descriptive properties of each traces
%   is then calculated (SNR, D/A correlation, etc).  Using defined
%   criteria, the user can then select a portion of the dataset.
%
%   Use "Batch Mode" to load all traces within a directory.
%   Be careful not to mix multiple experiments in the same directory!
%
%   Use "Save Traces" to save the resulting corrected+filtered traces.
%   This will also produce a log file with useful information.
%
%   The typical criteria used are:
%     FRET-lifetime    > 15 frames       (trace must show FRET)
%     D/A correlation  < 0.5             (remove aggregates)
%     Signal-to-Noise  > 8               (sufficient resolution)
%     Background noise < 1500            (background drift)
%     N. Donor Blinks  < 3               (remove aggregates)
%     Remove traces with multiple dyes = YES
%
%   NOTE that many of these criteria can bias the data, especially
%   correlation.  When comparing datasets, use the same criteria.

%   A gamma (sensitivity/quantum yield ratio) correction is used in
%   calculating total intensity and SNR.  The value comes from
%   cascadeConstants and was calculated for our equipment with ribosome
%   samples.  The value will vary based on equipment and sample studied.
%   
%   Fluorescence data are no longer stored in handles because loading
%   too many traces at once cause out-of-memory errors.  They are now
%   oaded on-demand through GetPickedTraces.  Post-synchronization of traces
%   is no longer implemented!
%   8/2007  -DT
%
%   FRET lifetime, N. Donor Blinks, signal overlap detection criteria all
%   added from original version by JBM.
%   4/2008  -DT


% Last Modified by GUIDE v2.5 03-Mar-2009 14:55:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @realtimeAnalysis_settings_OpeningFcn, ...
    'gui_OutputFcn',  @realtimeAnalysis_settings_OutputFcn, ...
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






%#########################################################################
%------------------------- INITIALIZATION (GUI) -------------------------%
%#########################################################################


%----------INITIALIZATION OF THE GUI----------%
% --- Executes just before realtimeAnalysis_settings is made visible.
function realtimeAnalysis_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to realtimeAnalysis_settings (see VARARGIN)


% Choose default command line output for realtimeAnalysis_settings
handles.output=hObject;


% Leave everything alone if the program is already running.
% This initialization proceedure will confuse the program state.
% if isfield(handles,'criteria'),
%     disp('Autotrace is already running!');
%     return;
% end

%---- PROGRAM CONSTANTS
constants = cascadeConstants();

handles.constants = constants;


%---- INITIAL VALUES FOR PICKING CRITERIA
% These can be changed while running the program.
criteria.minTotalIntesity=5000;
criteria.maxTotalIntesity=15000;
criteria.minTotalLifetime=20;
criteria.maxTotalLifetime=1485;
criteria.minCorrelation=-1.1;
criteria.maxCorrelation=0.5;
criteria.minFret=0.125;
criteria.minSNR=8;
criteria.maxBackground=1500;
criteria.maxDonorBlinks = 4;
criteria.minFretLifetime = 15;
criteria.overlap = 1;

handles.criteria = criteria;


% Initialize input fields with values defined above.
set(handles.MeanTotalIntensityLow,'String',num2str(handles.criteria.minTotalIntesity));
set(handles.MeanTotalIntensityHigh,'String',num2str(handles.criteria.maxTotalIntesity));
set(handles.fretSlopeThresh,'String',num2str(handles.criteria.minFret));
set(handles.FluorescenceLifetime,'String',num2str(handles.criteria.minTotalLifetime));
set(handles.FluorescenceLifetimeHigh,'String',num2str(handles.criteria.maxTotalLifetime));
set(handles.lowThresh,'String',num2str(handles.criteria.minCorrelation));
set(handles.highThresh,'String',num2str(handles.criteria.maxCorrelation));
set(handles.SignalNoiseThresh,'String',num2str(handles.criteria.minSNR));
set(handles.BackgroundNoiseThresh,'String',num2str(handles.criteria.maxBackground));
set(handles.editNCross,'String',num2str(handles.criteria.maxDonorBlinks));
set(handles.editAccLife,'String',num2str(handles.criteria.minFretLifetime));

% set(handles.FRETBinSize,'String',num2str(handles.contour_bin_size));


% % Reset CHECKBOXES to all being CHECKED by default
% stateval = 1;
% 
% set(handles.MeanTotalIntensityBox,'Value', stateval);
% % set(handles.MeanAcceptorIntensityBox,'Value',stateval);
% set(handles.FluorescenceLifetimeBox,'Value',stateval);
% set(handles.CorrelationCoefficientBox,'Value',stateval);
% set(handles.SignalNoiseBox,'Value',stateval);
% set(handles.BackgroundNoiseBox,'Value',stateval);
% set(handles.fretThreshold,'Value',stateval);
% set(handles.checkboxNCross,'Value',stateval);
% set(handles.checkboxAccLife,'Value',stateval);
% 
% 
% % Enable TEXTBOXES from input by default
% stateval = 'on';
% 
% set(handles.MeanTotalIntensityLow,'Enable',stateval);
% set(handles.MeanTotalIntensityHigh,'Enable',stateval);
% set(handles.fretSlopeThresh,'Enable',stateval);
% % set(handles.MeanAcceptorIntensity,'Enable',stateval);
% set(handles.FluorescenceLifetime,'Enable',stateval);
% set(handles.FluorescenceLifetimeHigh,'Enable',stateval);
% set(handles.lowThresh,'Enable',stateval);
% set(handles.highThresh,'Enable',stateval);
% set(handles.editNCross,'Enable',stateval);
% set(handles.editAccLife,'Enable',stateval);


% Update handles structure
guidata(hObject,handles);

% END FUNCTION realtimeAnalysis_settings_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = realtimeAnalysis_settings_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1}=handles.output;



% END FUNCTION realtimeAnalysis_settings_OutputFcn


function close_Callback(src,event)
% Set figure as hidden (not Visible) so its properties can be
% accessed without the window being in the way.
set(gcf,'Visible','off');






%#########################################################################
%------------------------- INPUT FIELD HANDLING -------------------------%
%#########################################################################



%----------CHECKBOX FOR MEAN TOTAL INTENSITY CRITERIA----------
% --- Executes on button press in MeanTotalIntensityBox.
function MeanTotalIntensityBox_Callback(hObject, eventdata, handles)
% If box is checked, get the values input by the user.
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minTotalIntesity=str2double(get(handles.MeanTotalIntensityLow,'String'));
    handles.criteria.maxTotalIntesity=str2double(get(handles.MeanTotalIntensityHigh,'String'));
    set(handles.MeanTotalIntensityLow,'Enable','on');
    set(handles.MeanTotalIntensityHigh,'Enable','on');
else
    % If the box is unchecked, use values which will nullify the criterium.
    handles.criteria.minTotalIntesity=[];
    handles.criteria.maxTotalIntesity=[];
    set(handles.MeanTotalIntensityLow,'Enable','off');
    set(handles.MeanTotalIntensityHigh,'Enable','off');
end
guidata(hObject,handles);


%----------CHECKBOX FOR FRET THRESHOLD----------
% --- Executes on button press in fretThreshold.
function fretThreshold_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minFret=str2double(get(handles.fretSlopeThresh,'String'));
    set(handles.fretSlopeThresh,'Enable','on');
else
    handles.criteria.minFret=[];
    set(handles.fretSlopeThresh,'Enable','off');
end
guidata(hObject,handles);


%----------CHECKBOX FOR FLUORESCENCE LIFETIME CRITERIA----------
% --- Executes on button press in FluorescenceLifetimeBox.
function FluorescenceLifetimeBox_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minTotalLifetime=str2double(get(handles.FluorescenceLifetime,'String'));
    handles.criteria.maxTotalLifetime=str2double(get(handles.FluorescenceLifetimeHigh,'String'));
    set(handles.FluorescenceLifetime,'Enable','on');
    set(handles.FluorescenceLifetimeHigh,'Enable','on');
else
    handles.criteria.minTotalLifetime=[];
    handles.criteria.maxTotalLifetime=[];
    set(handles.FluorescenceLifetime,'Enable','off');
    set(handles.FluorescenceLifetimeHigh,'Enable','off');
end
guidata(hObject,handles);


%----------CHECKBOX FOR CORRELATION COEFFICIENT CRITERIA----------
% --- Executes on button press in CorrelationCoefficientBox.
function CorrelationCoefficientBox_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minCorrelation=str2double(get(handles.lowThresh,'String'));
    handles.criteria.maxCorrelation=str2double(get(handles.highThresh,'String'));
    set(handles.lowThresh,'Enable','on');
    set(handles.highThresh,'Enable','on');
else
    handles.criteria.minCorrelation=[];
    handles.criteria.maxCorrelation=[];
    set(handles.lowThresh,'Enable','off');
    set(handles.highThresh,'Enable','off');
end
guidata(hObject,handles);


%----------CHECKBOX FOR SIGNAL-TO-NOISE CRITERIA----------
% --- Executes on button press in SignalNoiseBox.
function SignalNoiseBox_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minSNR=str2double(get(handles.SignalNoiseThresh,'String'));
    set(handles.SignalNoiseThresh,'Enable','on');
else
    handles.criteria.minSNR=[];
    set(handles.SignalNoiseThresh,'Enable','off');
end
guidata(hObject,handles);


%----------CHECKBOX FOR BACKGROUND NOISE CRITERIA----------
% --- Executes on button press in BackgroundNoiseBox.
function BackgroundNoiseBox_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.maxBackground=str2double(get(handles.BackgroundNoiseThresh,'String'));
    set(handles.BackgroundNoiseThresh,'Enable','on');
else
    handles.criteria.maxBackground=[];
    set(handles.BackgroundNoiseThresh,'Enable','off');
end
guidata(hObject,handles);


% --- Executes on button press in checkboxNCross.
function checkboxNCross_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkboxNCross

% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.maxDonorBlinks=str2double(get(handles.editNCross,'String'));
    set(handles.editNCross,'Enable','on');
else
    handles.criteria.maxDonorBlinks=[];
    set(handles.editNCross,'Enable','off');
end
guidata(hObject,handles);



% --- Executes on button press in checkboxAccLife.
function checkboxAccLife_Callback(hObject, eventdata, handles)
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minFretLifetime = str2double(get(handles.editAccLife,'String'));
    set(handles.editAccLife,'Enable','on');
else
    handles.criteria.minFretLifetime=[];
    set(handles.editAccLife,'Enable','off');
end
guidata(hObject,handles);


function MeanTotalIntensityLow_Callback(hObject, eventdata, handles)
handles.criteria.minTotalIntesity=str2double(get(hObject,'String'));
guidata(hObject,handles);

function MeanTotalIntensityHigh_Callback(hObject, eventdata, handles)
handles.criteria.maxTotalIntesity=str2double(get(hObject,'String'));
guidata(hObject,handles);

function fretSlopeThresh_Callback(hObject, eventdata, handles)
handles.criteria.minFret=str2double(get(hObject,'String'));
guidata(hObject,handles);

function FluorescenceLifetime_Callback(hObject, eventdata, handles)
handles.criteria.minTotalLifetime=str2double(get(hObject,'String'));
guidata(hObject,handles);

function FluorescenceLifetimeHigh_Callback(hObject, eventdata, handles)
handles.criteria.maxTotalLifetime=str2double(get(hObject,'String'));
guidata(hObject,handles);

function lowThresh_Callback(hObject, eventdata, handles)
handles.criteria.minCorrelation=str2double(get(hObject,'String'));
guidata(hObject,handles);

function highThresh_Callback(hObject, eventdata, handles)
handles.criteria.maxCorrelation=str2double(get(hObject,'String'));
guidata(hObject,handles);

function SignalNoiseThresh_Callback(hObject, eventdata, handles)
handles.criteria.minSNR=str2double(get(hObject,'String'));
guidata(hObject,handles);

function BackgroundNoiseThresh_Callback(hObject, eventdata, handles)
handles.criteria.maxBackground=str2double(get(hObject,'String'));
guidata(hObject,handles);

function FRETBinSize_Callback(hObject, eventdata, handles)
handles.criteria.contour_bin_size=str2double(get(hObject,'String'));
guidata(hObject,handles);

function editNCross_Callback(hObject, eventdata, handles)
handles.criteria.maxDonorBlinks = str2double(get(hObject,'String'));
guidata(hObject,handles);

function editAccLife_Callback(hObject, eventdata, handles)
handles.criteria.minFretLifetime = str2double(get(hObject,'String'));
guidata(hObject,handles);

function chkOverlap_Callback(hObject, eventdata, handles)
handles.criteria.overlap = get(hObject,'Value');
guidata(hObject,handles);



function txtGettracesThreshold_Callback(hObject, eventdata, handles)
handles.gettracesThreshold = str2num( get(hObject,'String') );
guidata(hObject,handles);

function txtGettracesOverlap_Callback(hObject, eventdata, handles)
handles.gettracesOverlap = str2num( get(hObject,'String') );
guidata(hObject,handles);

function chkGettracesAlignment_Callback(hObject, eventdata, handles)
handles.gettracesAlign = get(hObject,'Value');
guidata(hObject,handles);



% --- Executes on button press in btnOkay.
function btnOkay_Callback(hObject, eventdata, handles)
% hObject    handle to btnOkay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Inform the main application window that settings have changed
% and to redo analysis.
% set(gcf,'Visible','off');
% 
% handles.realtimeAnalysis_Notify;



function close_Callback(hObject, eventdata, handles)
return;

