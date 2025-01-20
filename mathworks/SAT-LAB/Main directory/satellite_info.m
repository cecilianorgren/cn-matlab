function varargout = satellite_info(varargin)
% SATELLITE_INFO MATLAB code for satellite_info.fig
%      SATELLITE_INFO, by itself, creates a new SATELLITE_INFO or raises the existing
%      singleton*.
%
%      H = SATELLITE_INFO returns the handle to a new SATELLITE_INFO or the handle to
%      the existing singleton*.
%
%      SATELLITE_INFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SATELLITE_INFO.M with the given input arguments.
%
%      SATELLITE_INFO('Property','Value',...) creates a new SATELLITE_INFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before satellite_info_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to satellite_info_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help satellite_info

% Last Modified by GUIDE v2.5 02-Feb-2017 13:44:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @satellite_info_OpeningFcn, ...
    'gui_OutputFcn',  @satellite_info_OutputFcn, ...
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


% --- Executes just before satellite_info is made visible.
function satellite_info_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to satellite_info (see VARARGIN)

% Choose default command line output for satellite_info
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes satellite_info wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%Declare global variables
global sat_data_ori sel_num

%Move GUI at the center of screen
movegui(gcf,'center')

if isempty(sat_data_ori) == 0
    
    sat_name = extractfield(sat_data_ori,'satellite_name');
    sat_name = sat_name(:);
    
    set(handles.listbox1,'String',sat_name)
    
    sel_num = 1;
    
    set(handles.edit1,'String',sat_data_ori(sel_num).satellite_name)
    set(handles.edit2,'String',sat_data_ori(sel_num).satellite_number_1)
    set(handles.edit3,'String',sat_data_ori(sel_num).classification)
    set(handles.edit4,'String',sat_data_ori(sel_num).launch_year)
    set(handles.edit5,'String',sat_data_ori(sel_num).launch_number)
    set(handles.edit6,'String',sat_data_ori(sel_num).launch_piece)
    set(handles.edit7,'String',sat_data_ori(sel_num).epoch_year)
    set(handles.edit8,'String',sat_data_ori(sel_num).epoch_day)
    set(handles.edit9,'String',sat_data_ori(sel_num).first_derivative_mean_motion)
    set(handles.edit10,'String',sat_data_ori(sel_num).second_derivative_mean_motion)
    set(handles.edit11,'String',sat_data_ori(sel_num).BSTAR_drag_term)
    set(handles.edit12,'String',sat_data_ori(sel_num).ephemeris_type)
    set(handles.edit13,'String',sat_data_ori(sel_num).element_number)
    set(handles.edit14,'String',sat_data_ori(sel_num).inclination)
    set(handles.edit15,'String',sat_data_ori(sel_num).right_ascension_AN)
    set(handles.edit16,'String',sat_data_ori(sel_num).eccentricity)
    set(handles.edit17,'String',sat_data_ori(sel_num).argument_of_perigee)
    set(handles.edit18,'String',sat_data_ori(sel_num).mean_anomaly)
    set(handles.edit19,'String',sat_data_ori(sel_num).mean_motion)
    set(handles.edit20,'String',sat_data_ori(sel_num).revolution_number)
    
end


% --- Outputs from this function are returned to the command line.
function varargout = satellite_info_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(~, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(~, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(~, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(~, ~, ~)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, ~, ~)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(~, ~, ~)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, ~, ~)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(~, ~, ~)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, ~, ~)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(~, ~, ~)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, ~, ~)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(~, ~, ~)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, ~, ~)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(~, ~, ~)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, ~, ~)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(~, ~, ~)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, ~, ~)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(~, ~, ~)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, ~, ~)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(~, ~, ~)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, ~, ~)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(~, ~, ~)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, ~, ~)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(~, ~, ~)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, ~, ~)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(~, ~, ~)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, ~, ~)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(~, ~, ~)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, ~, ~)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(~, ~, ~)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, ~, ~)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(~, ~, ~)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, ~, ~)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(~, ~, ~)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, ~, ~)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(~, ~, ~)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, ~, ~)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(~, ~, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
global sat_data_ori sat_data_new sel_num

if isempty(sat_data_ori) == 0

sel_num = get(handles.listbox1,'Value');

if get(handles.radiobutton1,'Value') == 1
    
    set(handles.edit1,'String',sat_data_ori(sel_num).satellite_name)
    set(handles.edit2,'String',sat_data_ori(sel_num).satellite_number_1)
    set(handles.edit3,'String',sat_data_ori(sel_num).classification)
    set(handles.edit4,'String',sat_data_ori(sel_num).launch_year)
    set(handles.edit5,'String',sat_data_ori(sel_num).launch_number)
    set(handles.edit6,'String',sat_data_ori(sel_num).launch_piece)
    set(handles.edit7,'String',sat_data_ori(sel_num).epoch_year)
    set(handles.edit8,'String',sat_data_ori(sel_num).epoch_day)
    set(handles.edit9,'String',sat_data_ori(sel_num).first_derivative_mean_motion)
    set(handles.edit10,'String',sat_data_ori(sel_num).second_derivative_mean_motion)
    set(handles.edit11,'String',sat_data_ori(sel_num).BSTAR_drag_term)
    set(handles.edit12,'String',sat_data_ori(sel_num).ephemeris_type)
    set(handles.edit13,'String',sat_data_ori(sel_num).element_number)
    set(handles.edit14,'String',sat_data_ori(sel_num).inclination)
    set(handles.edit15,'String',sat_data_ori(sel_num).right_ascension_AN)
    set(handles.edit16,'String',sat_data_ori(sel_num).eccentricity)
    set(handles.edit17,'String',sat_data_ori(sel_num).argument_of_perigee)
    set(handles.edit18,'String',sat_data_ori(sel_num).mean_anomaly)
    set(handles.edit19,'String',sat_data_ori(sel_num).mean_motion)
    set(handles.edit20,'String',sat_data_ori(sel_num).revolution_number)
    
elseif get(handles.radiobutton2,'Value') == 1
    
    set(handles.edit1,'String',sat_data_new(sel_num).satellite_name)
    set(handles.edit2,'String',sat_data_new(sel_num).satellite_number_1)
    set(handles.edit3,'String',sat_data_new(sel_num).classification)
    set(handles.edit4,'String',sat_data_new(sel_num).launch_year)
    set(handles.edit5,'String',sat_data_new(sel_num).launch_number)
    set(handles.edit6,'String',sat_data_new(sel_num).launch_piece)
    set(handles.edit7,'String',sat_data_new(sel_num).epoch_year)
    set(handles.edit8,'String',sat_data_new(sel_num).epoch_day)
    set(handles.edit9,'String',sat_data_new(sel_num).first_derivative_mean_motion)
    set(handles.edit10,'String',sat_data_new(sel_num).second_derivative_mean_motion)
    set(handles.edit11,'String',sat_data_new(sel_num).BSTAR_drag_term)
    set(handles.edit12,'String',sat_data_new(sel_num).ephemeris_type)
    set(handles.edit13,'String',sat_data_new(sel_num).element_number)
    set(handles.edit14,'String',sat_data_new(sel_num).inclination)
    set(handles.edit15,'String',sat_data_new(sel_num).right_ascension_AN)
    set(handles.edit16,'String',sat_data_new(sel_num).eccentricity)
    set(handles.edit17,'String',sat_data_new(sel_num).argument_of_perigee)
    set(handles.edit18,'String',sat_data_new(sel_num).mean_anomaly)
    set(handles.edit19,'String',sat_data_new(sel_num).mean_motion)
    set(handles.edit20,'String',sat_data_new(sel_num).revolution_number)
    
end

end


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, ~, ~)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(~, ~, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
global sat_data_ori sel_num

set(handles.radiobutton1,'Value',1)
set(handles.radiobutton2,'Value',0)

if isempty(sat_data_ori) == 0

set(handles.edit1,'String',sat_data_ori(sel_num).satellite_name)
set(handles.edit2,'String',sat_data_ori(sel_num).satellite_number_1)
set(handles.edit3,'String',sat_data_ori(sel_num).classification)
set(handles.edit4,'String',sat_data_ori(sel_num).launch_year)
set(handles.edit5,'String',sat_data_ori(sel_num).launch_number)
set(handles.edit6,'String',sat_data_ori(sel_num).launch_piece)
set(handles.edit7,'String',sat_data_ori(sel_num).epoch_year)
set(handles.edit8,'String',sat_data_ori(sel_num).epoch_day)
set(handles.edit9,'String',sat_data_ori(sel_num).first_derivative_mean_motion)
set(handles.edit10,'String',sat_data_ori(sel_num).second_derivative_mean_motion)
set(handles.edit11,'String',sat_data_ori(sel_num).BSTAR_drag_term)
set(handles.edit12,'String',sat_data_ori(sel_num).ephemeris_type)
set(handles.edit13,'String',sat_data_ori(sel_num).element_number)
set(handles.edit14,'String',sat_data_ori(sel_num).inclination)
set(handles.edit15,'String',sat_data_ori(sel_num).right_ascension_AN)
set(handles.edit16,'String',sat_data_ori(sel_num).eccentricity)
set(handles.edit17,'String',sat_data_ori(sel_num).argument_of_perigee)
set(handles.edit18,'String',sat_data_ori(sel_num).mean_anomaly)
set(handles.edit19,'String',sat_data_ori(sel_num).mean_motion)
set(handles.edit20,'String',sat_data_ori(sel_num).revolution_number)

end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(~, ~, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
global sat_data_ori sat_data_new sel_num

set(handles.radiobutton1,'Value',0)
set(handles.radiobutton2,'Value',1)

if isempty(sat_data_ori) == 0

set(handles.edit1,'String',sat_data_new(sel_num).satellite_name)
set(handles.edit2,'String',sat_data_new(sel_num).satellite_number_1)
set(handles.edit3,'String',sat_data_new(sel_num).classification)
set(handles.edit4,'String',sat_data_new(sel_num).launch_year)
set(handles.edit5,'String',sat_data_new(sel_num).launch_number)
set(handles.edit6,'String',sat_data_new(sel_num).launch_piece)
set(handles.edit7,'String',sat_data_new(sel_num).epoch_year)
set(handles.edit8,'String',sat_data_new(sel_num).epoch_day)
set(handles.edit9,'String',sat_data_new(sel_num).first_derivative_mean_motion)
set(handles.edit10,'String',sat_data_new(sel_num).second_derivative_mean_motion)
set(handles.edit11,'String',sat_data_new(sel_num).BSTAR_drag_term)
set(handles.edit12,'String',sat_data_new(sel_num).ephemeris_type)
set(handles.edit13,'String',sat_data_new(sel_num).element_number)
set(handles.edit14,'String',sat_data_new(sel_num).inclination)
set(handles.edit15,'String',sat_data_new(sel_num).right_ascension_AN)
set(handles.edit16,'String',sat_data_new(sel_num).eccentricity)
set(handles.edit17,'String',sat_data_new(sel_num).argument_of_perigee)
set(handles.edit18,'String',sat_data_new(sel_num).mean_anomaly)
set(handles.edit19,'String',sat_data_new(sel_num).mean_motion)
set(handles.edit20,'String',sat_data_new(sel_num).revolution_number)

end
