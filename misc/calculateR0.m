function varargout = calculateR0(varargin)
% CALCULATER0 MATLAB code for calculateR0.fig
%      CALCULATER0, by itself, creates a new CALCULATER0 or raises the existing
%      singleton*.
%
%      H = CALCULATER0 returns the handle to a new CALCULATER0 or the handle to
%      the existing singleton*.
%
%      CALCULATER0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALCULATER0.M with the given input arguments.
%
%      CALCULATER0('Property','Value',...) creates a new CALCULATER0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calculateR0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calculateR0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calculateR0

% Last Modified by GUIDE v2.5 11-Jan-2024 15:10:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calculateR0_OpeningFcn, ...
                   'gui_OutputFcn',  @calculateR0_OutputFcn, ...
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


% --- Executes just before calculateR0 is made visible.
function calculateR0_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calculateR0 (see VARARGIN)

% Choose default command line output for calculateR0
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes calculateR0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = calculateR0_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in btnBrowseDonor.
function btnBrowseDonor_Callback(~, ~, handles)
[f,p] = uigetfile( '*.txt','Acceptor Spectrum',get(handles.edDonorEmission,'String') );
if ~isequal(f,0)
    set( handles.edDonorEmission, 'String',fullfile(p,f) );
end



% --- Executes on button press in btnBrowseAcceptor.
function btnBrowseAcceptor_Callback(~, ~, handles)
[f,p] = uigetfile( '*.txt','Acceptor Spectrum',get(handles.edAcceptorAbsorbance,'String') );
if ~isequal(f,0)
    set( handles.edAcceptorAbsorbance, 'String',fullfile(p,f) );
end



% --- Executes on button press in btnCalcR0.
function btnCalcR0_Callback(~, ~, handles) %#ok<*DEFNU>
% Verify field values and calculate R0.

try
    % Load parameters from UI
    q   = str2double( get(handles.edDonorQY,'String') );
    ext = str2double( get(handles.edExtinctionCoef,'String') );
    n   = str2double( get(handles.edRefractiveIndex,'String') );
    K2  = str2double( get(handles.edKappa,'String') );
    
    wavelength = (400:1:900)';
    f = loadSpectrum( get(handles.edDonorEmission,'String'), wavelength );
    a = loadSpectrum( get(handles.edAcceptorAbsorbance,'String'), wavelength );
    
    plot( handles.axSpectra, wavelength, [f a] );
    ylim([0 1]);
    xlim([wavelength(1) wavelength(end)]);
    xlabel('Wavelength (nm)');

    auc = min(f,a);
    patch(wavelength,auc,'y');
    
    d_lambda = wavelength(2)-wavelength(1);
    assert( all( diff(wavelength)==d_lambda ) );

    % Calculate RO
    J = sum( (f/sum(f)) .* (ext*a) .* (wavelength.^4) );
    R0 = 0.211 * ( q * K2 * n^(-4) * J  )^(1/6);
    set( handles.edR0, 'String',num2str(R0,4) );

    %J = sum( (f/sum(f)) .* param.excoef.*(a/max(a)) .* (wavelength.^4) );
    %R0 = 0.211 * ( param.q * param.K2 * param.n^(-4) * J  )^(1/6);

    % FRET-distance relationship
    %subplot(1,2,2);
    %distance = 20:0.5:100;
    %fret = 1./( 1+ (distance/R0).^6 );
    %hplot = plot(distance, fret);
    %xlabel('Inter-dye distance (A)');
    %ylabel('FRET Efficiency');
    %datatip(hplot,R0,0.5);
    
catch e
    msgbox( ['Error: ' e.message] );
end

%end function btnCalcR0_Callback


% Load spectrum from file, interpolating to input wavelength values
function values = loadSpectrum(filename, wavelengths)

data = load(filename,'-ASCII');
assert( size(data,2)==2, 'Invalid file: must have two columns (wavelength and spectrum)' );
wavelengths = to_col(wavelengths);
values = interp1(data(:,1), data(:,2), wavelengths, 'linear',0);

% Crude baseline subtraction and normalization
values = values - min(values(values>0));
values = values / max(values);

%end function loadSpectrum



