function varargout = catalog_selection(varargin)
% CATALOG_SELECTION MATLAB code for catalog_selection.fig
%      CATALOG_SELECTION, by itself, creates a new CATALOG_SELECTION or raises the existing
%      singleton*.
%
%      H = CATALOG_SELECTION returns the handle to a new CATALOG_SELECTION or the handle to
%      the existing singleton*.
%
%      CATALOG_SELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CATALOG_SELECTION.M with the given input arguments.
%
%      CATALOG_SELECTION('Property','Value',...) creates a new CATALOG_SELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before catalog_selection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to catalog_selection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help catalog_selection

% Last Modified by GUIDE v2.5 02-Feb-2017 13:47:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @catalog_selection_OpeningFcn, ...
    'gui_OutputFcn',  @catalog_selection_OutputFcn, ...
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


% --- Executes just before catalog_selection is made visible.
function catalog_selection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to catalog_selection (see VARARGIN)

% Choose default command line output for catalog_selection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes catalog_selection wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%Move GUI at the center of screen
movegui(gcf,'center')


% --- Outputs from this function are returned to the command line.
function varargout = catalog_selection_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox35.
function checkbox35_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox35


% --- Executes on button press in checkbox36.
function checkbox36_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox36


% --- Executes on button press in checkbox37.
function checkbox37_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox37


% --- Executes on button press in checkbox38.
function checkbox38_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox38


% --- Executes on button press in checkbox39.
function checkbox39_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox39


% --- Executes on button press in checkbox40.
function checkbox40_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox40


% --- Executes on button press in checkbox41.
function checkbox41_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox41


% --- Executes on button press in checkbox31.
function checkbox31_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox31


% --- Executes on button press in checkbox32.
function checkbox32_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox32


% --- Executes on button press in checkbox33.
function checkbox33_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox33


% --- Executes on button press in checkbox34.
function checkbox34_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox34


% --- Executes on button press in checkbox27.
function checkbox27_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox27


% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox28


% --- Executes on button press in checkbox29.
function checkbox29_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox29


% --- Executes on button press in checkbox30.
function checkbox30_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox30


% --- Executes on button press in checkbox20.
function checkbox20_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox20


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21


% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22


% --- Executes on button press in checkbox23.
function checkbox23_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox23


% --- Executes on button press in checkbox24.
function checkbox24_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox24


% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25


% --- Executes on button press in checkbox26.
function checkbox26_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox26


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13


% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox14


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15


% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16


% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17


% --- Executes on button press in checkbox18.
function checkbox18_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox18


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global catalog sat_data_ori sat_data_new catalog_selection_check SAT_LAB_dir handles2

button_ans = questdlg('All files in Current Catalogs folder will be deleted. Do you want to continue?','Warning','Yes','No','Yes');

if strcmp(button_ans,'Yes') == 1
    
    catalog_selection_check = 0;
    
    catalog = zeros(41,2);
    
    %Make Current Catalogs folder
    warning('off','all')
    mkdir([SAT_LAB_dir 'Main directory\Current Catalogs'])
    warning('on','all')
    
    %Add Current Catalogs to path
    addpath([SAT_LAB_dir 'Main directory\Current Catalogs'])
    
    %Delete all files from Current Catalogs
    delete([SAT_LAB_dir 'Main directory\Current Catalogs\*'])
    
    %Add waitbar
    w_bar = waitbar(0,'Downloading selected Catalogs. Please wait...');
    
    %Download the Catalogs and update the waitbar
    
    if get(handles.checkbox1,'Value') == 1
        catalog(1,1) = 1;
        urlwrite('https://celestrak.com/NORAD/elements/weather.txt','weather.txt');
        movefile([cd '\weather.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\weather.txt']);
    end
    
    waitbar(1/41)
    
    if get(handles.checkbox2,'Value') == 1
        catalog(2,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/noaa.txt','noaa.txt');
        movefile([cd '\noaa.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\noaa.txt']);
    end
    waitbar(2/41)
    
    if get(handles.checkbox3,'Value') == 1
        catalog(3,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/goes.txt','goes.txt');
        movefile([cd '\goes.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\goes.txt']);
    end
    
    waitbar(3/41)
    
    if get(handles.checkbox4,'Value') == 1
        catalog(4,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/resource.txt','resource.txt');
        movefile([cd '\resource.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\resource.txt']);
    end
    
    waitbar(4/41)
    
    if get(handles.checkbox5,'Value') == 1
        catalog(5,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/sarsat.txt','sarsat.txt');
        movefile([cd '\sarsat.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\sarsat.txt']);
    end
    
    waitbar(5/41)
    
    if get(handles.checkbox6,'Value') == 1
        catalog(6,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/dmc.txt','dmc.txt');
        movefile([cd '\dmc.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\dmc.txt']);
    end
    
    waitbar(6/41)
    
    if get(handles.checkbox7,'Value') == 1
        catalog(7,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/tdrss.txt','tdrss.txt');
        movefile([cd '\tdrss.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\tdrss.txt']);
    end
    
    waitbar(7/41)
    
    if get(handles.checkbox8,'Value') == 1
        catalog(8,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/argos.txt','argos.txt');
        movefile([cd '\argos.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\argos.txt']);
    end
    
    waitbar(8/41)
    
    if get(handles.checkbox9,'Value') == 1
        catalog(9,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/geo.txt','geo.txt');
        movefile([cd '\geo.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\geo.txt']);
    end
    
    waitbar(9/41)
    
    if get(handles.checkbox10,'Value') == 1
        catalog(10,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/intelsat.txt','intelsat.txt');
        movefile([cd '\intelsat.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\intelsat.txt']);
    end
    
    waitbar(10/41)
    
    if get(handles.checkbox11,'Value') == 1
        catalog(11,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/gorizont.txt','gorizont.txt');
        movefile([cd '\gorizont.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\gorizont.txt']);
    end
    
    waitbar(11/41)
    
    if get(handles.checkbox12,'Value') == 1
        catalog(12,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/raduga.txt','raduga.txt');
        movefile([cd '\raduga.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\raduga.txt']);
    end
    
    waitbar(12/41)
    
    if get(handles.checkbox13,'Value') == 1
        catalog(13,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/molniya.txt','molniya.txt');
        movefile([cd '\molniya.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\molniya.txt']);
    end
    
    waitbar(13/41)
    
    if get(handles.checkbox14,'Value') == 1
        catalog(14,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/iridium.txt','iridium.txt');
        movefile([cd '\iridium.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\iridium.txt']);
    end
    
    waitbar(14/41)
    
    if get(handles.checkbox15,'Value') == 1
        catalog(15,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/orbcomm.txt','orbcomm.txt');
        movefile([cd '\orbcomm.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\orbcomm.txt']);
    end
    
    waitbar(15/41)
    
    if get(handles.checkbox16,'Value') == 1
        catalog(16,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/globalstar.txt','globalstar.txt');
        movefile([cd '\globalstar.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\globalstar.txt']);
    end
    
    waitbar(16/41)
    
    if get(handles.checkbox17,'Value') == 1
        catalog(17,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/amateur.txt','amateur.txt');
        movefile([cd '\amateur.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\amateur.txt']);
    end
    
    waitbar(17/41)
    
    if get(handles.checkbox18,'Value') == 1
        catalog(18,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/x-comm.txt','x-comm.txt');
        movefile([cd '\x-comm.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\x-comm.txt']);
    end
    
    waitbar(18/41)
    
    if get(handles.checkbox19,'Value') == 1
        catalog(19,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/other-comm.txt','other-comm.txt');
        movefile([cd '\other-comm.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\other-comm.txt']);
    end
    
    waitbar(19/41)
    
    if get(handles.checkbox20,'Value') == 1
        catalog(20,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/gps-ops.txt','gps-ops.txt');
        movefile([cd '\gps-ops.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\gps-ops.txt']);
    end
    
    waitbar(20/41)
    
    if get(handles.checkbox21,'Value') == 1
        catalog(21,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/glo-ops.txt','glo-ops.txt');
        movefile([cd '\glo-ops.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\glo-ops.txt']);
    end
    
    waitbar(21/41)
    
    if get(handles.checkbox22,'Value') == 1
        catalog(22,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/galileo.txt','galileo.txt');
        movefile([cd '\galileo.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\galileo.txt']);
    end
    
    waitbar(22/41)
    
    if get(handles.checkbox23,'Value') == 1
        catalog(23,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/beidou.txt','beidou.txt');
        movefile([cd '\beidou.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\beidou.txt']);
    end
    
    waitbar(23/41)
    
    if get(handles.checkbox24,'Value') == 1
        catalog(24,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/sbas.txt','sbas.txt');
        movefile([cd '\sbas.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\sbas.txt']);
    end
    
    waitbar(24/41)
    
    if get(handles.checkbox25,'Value') == 1
        catalog(25,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/nnss.txt','nnss.txt');
        movefile([cd '\nnss.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\nnss.txt']);
    end
    
    waitbar(25/41)
    
    if get(handles.checkbox26,'Value') == 1
        catalog(26,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/musson.txt','musson.txt');
        movefile([cd '\musson.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\musson.txt']);
    end
    
    waitbar(26/41)
    
    if get(handles.checkbox27,'Value') == 1
        catalog(27,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/science.txt','science.txt');
        movefile([cd '\science.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\science.txt']);
    end
    
    waitbar(27/41)
    
    if get(handles.checkbox28,'Value') == 1
        catalog(28,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/geodetic.txt','geodetic.txt');
        movefile([cd '\geodetic.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\geodetic.txt']);
    end
    
    waitbar(28/41)
    
    if get(handles.checkbox29,'Value') == 1
        catalog(29,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/engineering.txt','engineering.txt');
        movefile([cd '\engineering.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\engineering.txt']);
    end
    
    waitbar(29/41)
    
    if get(handles.checkbox30,'Value') == 1
        catalog(30,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/education.txt','education.txt');
        movefile([cd '\education.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\education.txt']);
    end
    
    waitbar(30/41)
    
    if get(handles.checkbox31,'Value') == 1
        catalog(31,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/military.txt','military.txt');
        movefile([cd '\military.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\military.txt']);
    end
    
    waitbar(31/41)
    
    if get(handles.checkbox32,'Value') == 1
        catalog(32,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/radar.txt','radar.txt');
        movefile([cd '\radar.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\radar.txt']);
    end
    
    waitbar(32/41)
    
    if get(handles.checkbox33,'Value') == 1
        catalog(33,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/cubesat.txt','cubesat.txt');
        movefile([cd '\cubesat.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\cubesat.txt']);
    end
    
    waitbar(33/41)
    
    if get(handles.checkbox34,'Value') == 1
        catalog(34,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/other.txt','other.txt');
        movefile([cd '\other.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\other.txt']);
    end
    
    waitbar(34/41)
    
    if get(handles.checkbox35,'Value') == 1
        catalog(35,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/tle-new.txt','tle-new.txt');
        movefile([cd '\tle-new.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\tle-new.txt']);
    end
    
    waitbar(35/41)
    
    if get(handles.checkbox36,'Value') == 1
        catalog(36,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/stations.txt','stations.txt');
        movefile([cd '\stations.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\stations.txt']);
    end
    
    waitbar(36/41)
    
    if get(handles.checkbox37,'Value') == 1
        catalog(37,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/visual.txt','visual.txt');
        movefile([cd '\visual.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\visual.txt']);
    end
    
    waitbar(37/41)
    
    if get(handles.checkbox38,'Value') == 1
        catalog(38,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/1999-025.txt','1999-025.txt');
        movefile([cd '\1999-025.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\1999-025.txt']);
    end
    
    waitbar(38/41)
    
    if get(handles.checkbox39,'Value') == 1
        catalog(39,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/iridium-33-debris.txt','iridium-33-debris.txt');
        movefile([cd '\iridium-33-debris.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\iridium-33-debris.txt']);
    end
    
    waitbar(39/41)
    
    if get(handles.checkbox40,'Value') == 1
        catalog(40,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/cosmos-2251-debris.txt','cosmos-2251-debris.txt');
        movefile([cd '\cosmos-2251-debris.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\cosmos-2251-debris.txt']);
    end
    
    waitbar(40/41)
    
    if get(handles.checkbox41,'Value') == 1
        catalog(41,1) = 1;
        urlwrite('http://celestrak.com/NORAD/elements/2012-044.txt','2012-044.txt');
        movefile([cd '\2012-044.txt'],[SAT_LAB_dir 'Main directory\Current Catalogs\2012-044.txt']);
    end
    
    waitbar(41/41)
    
    %Close the waitbar
    close(w_bar)
    
    if max(max(catalog)) == 0
        
        warndlg('Please select at least one satellite catalog.')
        
    else
        
        %Copy downloaded files to Saved Catalogs
        copyfile([SAT_LAB_dir 'Main directory\Current Catalogs\*.txt'],[SAT_LAB_dir 'Main directory\Saved Catalogs\']);
        
        %Get all Catalog filenames
        file_names = dir([SAT_LAB_dir 'Main directory\Current Catalogs\*.txt']);
        
        %Add a waitbar
        w_bar = waitbar(0,'Reading selected Catalogs. Please wait...');
        
        %Read all Catalog files and update the waitbar
        sat_num = zeros(size(file_names));
        
        for i = 1:size(file_names,1);
            
            current_file = file_names(i,1).name;
            
            [cur_data_ori,cur_data_new] = f_read_tle(current_file);
            
            %Get number of satellites of current TLE file
            sat_num(i,1) = size(cur_data_ori,1);
            
            if i == 1
                
                sat_data_ori = cur_data_ori;
                sat_data_new = cur_data_new;
            else
                
                sat_data_ori = [sat_data_ori;cur_data_ori];
                sat_data_new = [sat_data_new;cur_data_new];
                
            end
            
            waitbar(i/size(file_names,1))
            
        end
        
        idx = catalog(:,1) > 0;
        
        catalog(idx,2) = sat_num;
        
        %Close the waitbar
        close(w_bar)
        
        %Close form
        close catalog_selection
        
        catalog_selection_check = 1;
        
        %Import names in the listbox
        sat_name = extractfield(sat_data_ori,'satellite_name');
        sat_name = sat_name(:);
        
        satlab_h   = findobj('Tag','satlab');
        satlab_gui = guidata(satlab_h);
        
        def_str = get(satlab_gui.popupmenu1,'String');
        
        if iscell(def_str) == 1
            
            def_name = def_str{1,1};
            
        elseif iscell(def_str) == 0
            
            def_name = def_str;
            
        end
        
        set(satlab_gui.popupmenu1,'String',[def_name;sat_name])
        set(satlab_gui.popupmenu1,'Value',1)
        set(satlab_gui.checkbox11,'Value',0)
        satlab('checkbox11_Callback',hObject, eventdata, handles2)
        
    end
    
end
