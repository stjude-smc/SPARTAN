function varargout = realtimeAnalysis_settings(varargin)
% REALTIMEANALYSIS_SETTINGS  Trace processing and filtering
%


% Last Modified by GUIDE v2.5 31-Mar-2010 18:26:23

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

%warning off MATLAB:divideByZero


% Setup callback to calling object. When the user clicks "Save", the
% criteria are sent back through this callback.
if numel(varargin) > 0
    if ischar(varargin{1}), error('unknown argument: %s',varargin{1}); end
    handles.updateCallback = varargin{1};
    handles.updateObject   = varargin{2};
end


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
criteria = constants.defaultAutotraceCriteria;


%---- Setup drop-down boxes listing trace statistics.

% Get names of trace statistics
ln = traceStat;  %get long statistic names
handles.statLongNames = ln;
longNames  = struct2cell(ln);
shortNames = fieldnames(ln);

% Setup special criteria selection drop-down boxes.
handles.nCriteriaBoxes = 8;
criteriaNames = [{''}; longNames];

for id=1:handles.nCriteriaBoxes
    set( handles.(['cboCriteria' num2str(id)]), 'String', criteriaNames );
end

% Setup default values (if no criteria seen) for special boxes.
set( handles.chkOverlap, 'Value',0 );


%---- Set default selections for the drop-down boxes.
fnames = fieldnames(criteria);
i = 1; %index into GUI elements

for id=1:numel(fnames)
    % There is a seperate control for removing overlapped/contamined
    % traces. BUT eq_overlap=1 cannot be handled here and will be placed in
    % the dropdown boxes instead.
    if strcmp(fnames{id},'eq_overlap') && criteria.eq_overlap==0,
        set( handles.chkOverlap, 'Value',1 );
        continue;
    end
    
    % Determine the short name of the criteria
    sepIndx = strfind(fnames{id},'_');
    criteriaName = fnames{id}(sepIndx+1:end);
    r = strcmp( shortNames, criteriaName );
    criteriaIndex = find(r,1);
    assert( ~isempty(criteriaIndex), 'Invalid criteria name!' );
    
    % Determine the equality sign used.
    equalityName = fnames{id}(1:sepIndx-1);
    r = strcmp( {'min','max','eq'}, equalityName );
    equalityIndex = find( r, 1 );
    assert( ~isempty(equalityIndex), 'Invalid equality name!' );
    
    % Set the dropdown boxes and criteria value textbox
    set( handles.(['cboCriteria' num2str(i)]), 'Value',  criteriaIndex+1 );
    set( handles.(['cboEquality' num2str(i)]), 'Value',  equalityIndex+1 );
    set( handles.(['edCriteria'  num2str(i)]), 'String', ...
                                           num2str(criteria.(fnames{id})) );
                                       
    i = i+1;
end

% This is set at the end since it has a non-standard name...
handles.criteria = criteria;


%---- Initial values for gettraces options
options = constants.gettracesDefaultParams;
handles.options = options;

% Set field values in GUI
set( handles.txtGettracesThreshold,   'String',num2str(options.don_thresh)   );
set( handles.txtGettracesOverlap,     'String',num2str(options.overlap_thresh)     );
set( handles.txtGettracesIntegration, 'String',num2str(options.nPixelsToSum) );


%---- Update handles structure
guidata(hObject,handles);

% END FUNCTION autotrace_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = realtimeAnalysis_settings_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1}=handles.output;



% END FUNCTION realtimeAnalysis_settings_OutputFcn


function close_Callback(hObject)
% Set figure as hidden (not Visible) so its properties can be
% accessed without the window being in the way. The window is only destroyed
% when the realtimeAnalysis GUI window is closed.
set(hObject,'Visible','off');






%#########################################################################
%------------------------- INPUT FIELD HANDLING -------------------------%
%#########################################################################


% --- CALLED when any of the gettraces options textboxes are changed.
% Updates the settings structure.
function updateOptions_Callback( hObject, handles, optionName )

val = [];
if strcmpi( get(hObject,'Style'), 'checkbox' )
    val = get(hObject,'Value');
else
    if ~isempty( get(hObject,'String') ),
        val = str2double( get(hObject,'String') );
    end
end

handles.options.(optionName) = val;
guidata(hObject,handles);



% --- CALLED when any of the primary selection criteria checkbox setting is
% changed. Updates the selection criteria automatically. 
function criteriaCheckbox_Callback(hObject, handles, criteriaName, textboxName)

isChecked = (get(hObject,'Value')==get(hObject,'Max'));

% Set textbox state, if there is one associated with this checkbox.
if nargin>=4, %checkbox with associated textbox
    textbox = handles.(textboxName);

    states = {'off','on'};
    set(textbox,'Enable',states{isChecked+1});
    
    val = str2double(get(textbox,'String'));
    
else %just a checkbox
    val = isChecked;
end

% overlap is alittle different than the others. If the box is checked, we
% want the criteria to be "eq_overlap==0", but the "value" is 1.
if strcmp(criteriaName,'eq_overlap'),
    val = ~val;
end

% Add/remove criteria.
if isChecked
    handles.criteria(1).(criteriaName) = val;
elseif isfield(handles.criteria,criteriaName),
    handles.criteria = rmfield(handles.criteria,criteriaName);
end


guidata(hObject,handles);



%----------APPLIES PICKING CRITERIA TO TRACES----------
% --- Executes on button press in PickTraces.
function PickTraces_Callback(hObject, handles)
% Generate a selection criteria structure (see pickTraces.m) using the
% "Specialized Selection Criteria" drop-down boxes.

criteria = struct([]);

% Update criteria for combo-box selections
shortNames = fieldnames(handles.statLongNames);
equalityText = {'min_','max_','eq_'};

for id=1:handles.nCriteriaBoxes
    selection = get( handles.(['cboCriteria' num2str(id)]), 'Value' );
    if selection==1, continue; end %no selection
    
    equality = get(handles.(['cboEquality' num2str(id)]),'Value');
    if equality==1, continue; end %no inequality selected
    
    criteriaName = [equalityText{equality-1} shortNames{selection-1}];
    edText = get(handles.(['edCriteria' num2str(id)]),'String');
    if isempty(edText), continue; end %no criteria value.
    
    criteria(1).(criteriaName) = ...
        str2double( get(handles.(['edCriteria' num2str(id)]),'String') );
    
end

% Save criteria.
handles.criteria = criteria;
guidata(hObject,handles);

% Handle other criteria with checkboxes.
criteriaCheckbox_Callback( handles.chkOverlap, handles, 'eq_overlap' );
handles = guidata(hObject);
criteriaCheckbox_Callback( handles.chkTotalSigma, handles, 'maxTotalSigma', 'edIntSigma' );


% END FUNCTION PickTraces_Callback



% --- Executes on button press in btnOkay.
function btnOkay_Callback(hObject, eventdata, handles)
% hObject    handle to btnOkay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Inform the main application window that settings have changed
% and to redo analysis.

result.gettracesOptions = handles.options;
result.criteria = handles.criteria;

% disp( handles.options );
% disp( handles.criteria );

handles.updateCallback( result, handles.updateObject );


