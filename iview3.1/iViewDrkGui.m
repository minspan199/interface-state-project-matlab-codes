function varargout = iViewDrkGui(varargin)
% IVIEWDRKGUI M-file for iViewDrkGui.fig
%      IVIEWDRKGUI, by itself, creates a new IVIEWDRKGUI or raises the existing
%      singleton*.
%
%      H = IVIEWDRKGUI returns the handle to a new IVIEWDRKGUI or the handle to
%      the existing singleton*.
%
%      IVIEWDRKGUI('Property','Value',...) creates a new IVIEWDRKGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to iViewDrkGui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      IVIEWDRKGUI('CALLBACK') and IVIEWDRKGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in IVIEWDRKGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iViewDrkGui

% Last Modified by GUIDE v2.5 21-Apr-2005 14:25:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iViewDrkGui_OpeningFcn, ...
                   'gui_OutputFcn',  @iViewDrkGui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before iViewDrkGui is made visible.
function iViewDrkGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for iViewDrkGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes iViewDrkGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = iViewDrkGui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function threshDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in isLinAll.
function isLinAll_Callback(hObject, eventdata, handles)
% hObject    handle to isLinAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isLinAll


% --- Executes on button press in isLogAll.
function isLogAll_Callback(hObject, eventdata, handles)
% hObject    handle to isLogAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isLogAll


