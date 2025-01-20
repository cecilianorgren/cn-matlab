function varargout = satlab(varargin)
% SATLAB MATLAB code for satlab.fig
%      SATLAB, by itself, creates a new SATLAB or raises the existing
%      singleton*.
%
%      H = SATLAB returns the handle to a new SATLAB or the handle to
%      the existing singleton*.
%
%      SATLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SATLAB.M with the given input arguments.
%
%      SATLAB('Property','Value',...) creates a new SATLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before satlab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to satlab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help satlab

% Last Modified by GUIDE v2.5 30-Mar-2017 14:54:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @satlab_OpeningFcn, ...
    'gui_OutputFcn',  @satlab_OutputFcn, ...
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


% --- Executes just before satlab is made visible.
function satlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to satlab (see VARARGIN)

% Choose default command line output for satlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes satlab wait for user response (see UIRESUME)
% uiwait(handles.satlab);

%Clear global variables
clearvars -global -except topo_res coast_res day_night_res

%Declare global variables
global POS_EARTH POS_COAST X_EARTH Y_EARTH Z_EARTH  X_COAST Y_COAST Z_COAST POS_t LON LAT TOPOGRAPHY row col row2 col2 long_corr lat_corr irf_axes efrf_axes
global GM period_Earth omega_Earth O a e w W i M theta_z T_num T t era
global h1a h2a h3a h4a h5a h6a h7a h8a h9a h1b h2b h3b h4b h5b h6b h7b h8b h9b h1c h2c h3c h4c h9c
global earth_color coastline_color satellite_color orbit_color velocity_color radial_color irf_color efrf_color
global sat_data_ori sat_data_new real_time_timer current_time_timer LONGITUDE_t LATITUDE_t X_t Y_t Z_t day_night_map
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat phi_sat phi_sat_corr lamda_sat lamda_sat_corr
global SAT_LAB_dir

%Move GUI at the center of screen
movegui(gcf,'center')

%Set current directory
SAT_LAB_dir                                 = which('satlab.m');
SAT_LAB_dir                                 = SAT_LAB_dir(1:end-8);

addpath(genpath(SAT_LAB_dir))

%Load Earth points and coastline data
if get(handles.popupmenu2,'Value') == 1
    
    load EARTH_2
    
elseif get(handles.popupmenu2,'Value') == 2
    
    load EARTH_1
    
elseif get(handles.popupmenu2,'Value') == 3
    
    load EARTH_05
    
end

if get(handles.popupmenu3,'Value') == 1
    
    load COAST_3D
    
elseif get(handles.popupmenu3,'Value') == 2
    
    load COAST_C_3D
    
elseif get(handles.popupmenu3,'Value') == 3
    
    load COAST_L_3D
    
end

load('coast.mat')

%Initialize editboxes
set(handles.edit1,'String',num2str(25000000))
set(handles.edit2,'String',num2str(0.60))
set(handles.edit3,'String',num2str(10))
set(handles.edit4,'String',num2str(195))
set(handles.edit5,'String',num2str(52))
set(handles.edit6,'String',num2str(10))
set(handles.edit7,'String',num2str(10))
set(handles.edit8,'String',num2str(1))
set(handles.edit9,'String',num2str(20000000))
set(handles.edit10,'String','')
set(handles.edit11,'String','')
set(handles.edit12,'String','')

set(handles.edit10,'Enable','off')
set(handles.edit11,'Enable','off')
set(handles.edit12,'Enable','off')

%Initialize slider position
set(handles.slider1,'Value',str2double(get(handles.edit1,'String')))
set(handles.slider2,'Value',str2double(get(handles.edit2,'String')))
set(handles.slider3,'Value',str2double(get(handles.edit3,'String')))
set(handles.slider4,'Value',str2double(get(handles.edit4,'String')))
set(handles.slider5,'Value',str2double(get(handles.edit5,'String')))
set(handles.slider6,'Value',str2double(get(handles.edit6,'String')))
set(handles.slider7,'Value',str2double(get(handles.edit7,'String')))
set(handles.slider8,'Value',str2double(get(handles.edit8,'String')))
set(handles.slider9,'Value',str2double(get(handles.edit9,'String')))

%Initialize checkboxes
set(handles.checkbox1,'Value',1)
set(handles.checkbox2,'Value',1)
set(handles.checkbox3,'Value',1)
set(handles.checkbox4,'Value',1)
set(handles.checkbox5,'Value',1)
set(handles.checkbox6,'Value',1)
set(handles.checkbox7,'Value',1)
set(handles.checkbox8,'Value',1)
set(handles.checkbox9,'Value',0)
set(handles.checkbox10,'Value',0)
set(handles.checkbox11,'Value',0)
set(handles.checkbox12,'Value',0)

%Initialize popupmenu
set(handles.popupmenu1,'Value',1)

%Initialize other GUI properties
set(handles.text14,'String','')

if get(handles.checkbox1,'Value') == 0
    
    set(handles.edit9,'Enable','On')
    set(handles.slider9,'Enable','On')
    
elseif get(handles.checkbox1,'Value') == 1
    
    set(handles.edit9,'Enable','Off')
    set(handles.slider9,'Enable','Off')
    
end

if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
    
    set(handles.checkbox3,'Enable','on')
    
elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
    
    set(handles.checkbox3,'Enable','on')
    
elseif get(handles.checkbox2,'Value') == 0
    
    set(handles.checkbox3,'Enable','off')
    
end

if get(handles.checkbox11,'Value') == 0
    
    set(handles.checkbox12,'Value',0)
    set(handles.checkbox12,'Enable','off')
    
elseif get(handles.checkbox11,'Value') == 1
    
    set(handles.checkbox12,'Enable','on')
    
end

%Initialize colors
earth_color                                 = get(handles.pushbutton1,'BackGroundColor');
coastline_color                             = get(handles.pushbutton2,'BackGroundColor');
satellite_color                             = get(handles.pushbutton3,'BackGroundColor');
orbit_color                                 = get(handles.pushbutton4,'BackGroundColor');
velocity_color                              = get(handles.pushbutton5,'BackGroundColor');
radial_color                                = get(handles.pushbutton6,'BackGroundColor');
irf_color                                   = get(handles.pushbutton7,'BackGroundColor');
efrf_color                                  = get(handles.pushbutton8,'BackGroundColor');

%Get all Catalog filenames
file_names                                  = dir([SAT_LAB_dir 'Main directory\Saved Catalogs\*.txt']);

%Add a waitbar
w_bar                                       = waitbar(0,'Reading saved Catalogs. Please wait...');

%Read all Catalog files
sat_num                                     = zeros(size(file_names));

for ii = 1:size(file_names,1);
    
    current_file                            = file_names(ii,1).name;
    [cur_data_ori,cur_data_new]             = f_read_tle(current_file);
    
    %Get number of satellites of current TLE file
    sat_num(ii,1)                           = size(cur_data_ori,1);
    
    if ii == 1
        
        sat_data_ori                        = cur_data_ori;
        sat_data_new                        = cur_data_new;
        
    else
        
        sat_data_ori                        = [sat_data_ori;cur_data_ori];
        sat_data_new                        = [sat_data_new;cur_data_new];
        
    end
    
    %Update the waitbar
    waitbar(ii/size(file_names,1))
    
end

%Close the waitbar
close(w_bar)

%Check is there are no Catalog files
if isempty(sat_data_ori) == 0
    
    sat_name                                = extractfield(sat_data_ori,'satellite_name');
    sat_name                                = sat_name(:);
    def_str                                 = get(handles.popupmenu1,'String');
    
    if iscell(def_str) == 1
        
        def_name                            = def_str{1,1};
        
    elseif iscell(def_str) == 0
        
        def_name                            = def_str;
        
    end
    
    set(handles.popupmenu1,'String',[def_name;sat_name])
    set(handles.popupmenu1,'Value',1)
    
end

%Set timers
real_time_timer                             = timer;
set(real_time_timer, 'ExecutionMode', 'fixedrate');
set(real_time_timer, 'Period', 5);
set(real_time_timer, 'TimerFcn', {@orbit_real_time, handles});

current_time_timer                          = timer;
set(current_time_timer, 'ExecutionMode', 'fixedrate');
set(current_time_timer, 'Period', 5);
set(current_time_timer, 'TimerFcn', {@show_current_time, handles});

start(current_time_timer)

%Set constants
GM                                          = 3.98199e+14;

%Initialize day night map
if get(handles.popupmenu4,'Value') == 1
    
    [LONGITUDE_t,LATITUDE_t]                = meshgrid(0:2:360,-90:2:90);
    
elseif get(handles.popupmenu4,'Value') == 2
    
    [LONGITUDE_t,LATITUDE_t]                = meshgrid(0:1:360,-90:1:90);
    
elseif get(handles.popupmenu4,'Value') == 3
    
    [LONGITUDE_t,LATITUDE_t]                = meshgrid(0:0.5:360,-90:0.5:90);
    
end

[X_t,Y_t,Z_t]                               = f_geo2cart(LATITUDE_t,LONGITUDE_t,10000,'GRS 1980');
[row2,col2]                                 = size(X_t);

x_t                                         = X_t(:);
y_t                                         = Y_t(:);
z_t                                         = Z_t(:);

POS_t                                       = [x_t';y_t';z_t'];
day_night_map                               = NaN(size(LONGITUDE_t));

%Calculate Earth's angular velocity
period_Earth                                = (23 + 56/60 + 4.0910/3600)*3600; %sec
omega_Earth                                 = 2*pi/period_Earth; %rad/sec

O                                           = [0             , -omega_Earth  , 0 ;
                                               omega_Earth   , 0             , 0 ;
                                               0             , 0             , 0];

%Process the Earth points and coastline data
[row,col]                                   = size(X_EARTH);

x_earth                                     = X_EARTH(:);
y_earth                                     = Y_EARTH(:);
z_earth                                     = Z_EARTH(:);

POS_EARTH                                   = [x_earth';y_earth';z_earth'];
POS_COAST                                   = [X_COAST';Y_COAST';Z_COAST'];

long(long < 0)                              = long(long < 0) + 360;

%Fix east-west lines bug
[lat_corr,long_corr]                        = f_fix_lines(lat,long,200);

%Propagate the orbit
a                                           = str2double(get(handles.edit1,'String'));
e                                           = str2double(get(handles.edit2,'String'));
w                                           = str2double(get(handles.edit3,'String'));
W                                           = str2double(get(handles.edit4,'String'));
i                                           = str2double(get(handles.edit5,'String'));
M                                           = str2double(get(handles.edit6,'String'));
theta_z                                     = str2double(get(handles.edit7,'String'))*pi/180;
T_num                                       = str2double(get(handles.edit8,'String'));
T                                           = 2*pi*sqrt(a^3/GM)/period_Earth;
t                                           = (0:T/360:T_num*T+T/360)';
era                                         = theta_z + omega_Earth*t*period_Earth;

%Inertial Reference Frame
[x_sat,y_sat,z_sat,v_x_sat,v_y_sat,v_z_sat] = f_kepl2svec(GM,a,e,w,W,i,M,0,t);

irf_axes                                    = [20000000,        0,        0 ;
                                                      0, 20000000,        0 ;
                                                      0,        0, 20000000];

%Initialize variables
X_sat                                       = zeros(size(x_sat));
Y_sat                                       = zeros(size(y_sat));
Z_sat                                       = zeros(size(z_sat));

V_X_sat                                     = zeros(size(v_x_sat));
V_Y_sat                                     = zeros(size(v_y_sat));
V_Z_sat                                     = zeros(size(v_z_sat));

%Earth Fixed Reference Frame
efrf_axes                                   = [20000000,        0,        0 ;
                                                      0, 20000000,        0 ;
                                                      0,        0, 20000000];

for ii = 1:size(x_sat,1)
    
    pos_vec                                 = [x_sat(ii,1);y_sat(ii,1);z_sat(ii,1)];
    vel_vec                                 = [v_x_sat(ii,1);v_y_sat(ii,1);v_z_sat(ii,1)];
    
    POS_VEC                                 = f_rotz(era(ii,1))*pos_vec;
    VEL_VEC                                 = f_rotz(era(ii,1))*vel_vec - O*POS_VEC;
    
    X_sat(ii,1)                             = POS_VEC(1,1);
    Y_sat(ii,1)                             = POS_VEC(2,1);
    Z_sat(ii,1)                             = POS_VEC(3,1);
    
    V_X_sat(ii,1)                           = VEL_VEC(1,1);
    V_Y_sat(ii,1)                           = VEL_VEC(2,1);
    V_Z_sat(ii,1)                           = VEL_VEC(3,1);
    
end

POS_EARTH_NEW                               = f_rotz(-theta_z)*POS_EARTH;
POS_COAST_NEW                               = f_rotz(-theta_z)*POS_COAST;

X_EARTH_NEW                                 = reshape(POS_EARTH_NEW(1,:),row,col);
Y_EARTH_NEW                                 = reshape(POS_EARTH_NEW(2,:),row,col);
Z_EARTH_NEW                                 = reshape(POS_EARTH_NEW(3,:),row,col);

X_COAST_NEW                                 = POS_COAST_NEW(1,:)';
Y_COAST_NEW                                 = POS_COAST_NEW(2,:)';
Z_COAST_NEW                                 = POS_COAST_NEW(3,:)';

IRF_AXES                                    = f_rotz(-theta_z)*irf_axes;
EFRF_AXES                                   = f_rotz(theta_z)*efrf_axes;

%Earth Fixed Reference Frame (Map Projection)
[phi_sat,lamda_sat,~]                       = f_cart2geo(X_sat,Y_sat,Z_sat,'GRS 1980');

%Fix east-west lines bug
[phi_sat_corr,lamda_sat_corr]               = f_fix_lines(phi_sat,lamda_sat,200);

%Plot to Inertial Reference Frame
axes(handles.axes1)

%Plot Earth
h1a                                         = surf(X_EARTH_NEW,Y_EARTH_NEW,Z_EARTH_NEW,TOPOGRAPHY);

if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
    
    demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
    shading interp
    
elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
    
    set(h1a,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color);
    
elseif get(handles.checkbox2,'Value') == 0
    
    set(h1a,'Visible','off')
    
end

hold on

%Plot day night map
h9a                                         = surf(X_t,Y_t,Z_t,day_night_map);

shading interp
set(h9a,'Visible','off')
alpha(h9a,0.5)

%Plot coastline
h2a                                         = plot3(X_COAST_NEW,Y_COAST_NEW,Z_COAST_NEW,'Color',coastline_color);

if get(handles.checkbox4,'Value') == 0
    
    set(h2a,'Visible','off')
    
end

%Plot satellite orbit
h3a                                         = plot3(x_sat(1:362),y_sat(1:362),z_sat(1:362),'Color',orbit_color,'Linewidth',2);

if get(handles.checkbox6,'Value') == 0
    
    set(h3a,'Visible','off')
    
end

%Plot satellite
h4a                                         = scatter3(x_sat(1,1),y_sat(1,1),z_sat(1,1),100,'filled','MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);

if get(handles.checkbox5,'Value') == 0
    
    set(h4a,'Visible','off')
    
end

%Plot satellite velocity
h5a                                         = quiver3(x_sat(1,1),y_sat(1,1),z_sat(1,1),v_x_sat(1,1)*1000,v_y_sat(1,1)*1000,v_z_sat(1,1)*1000,'Color',velocity_color,'Linewidth',2);

if get(handles.checkbox7,'Value') == 0
    
    set(h5a,'Visible','off')
    
end

%Plot satellite radial distance
h6a                                         = plot3([0 x_sat(1,1)],[0 y_sat(1,1)],[0 z_sat(1,1)],'Color',radial_color);

if get(handles.checkbox8,'Value') == 0
    
    set(h6a,'Visible','off')
    
end

%Plot Inertial Reference Frame axes
h7a                                         = quiver3([0;0;0],[0;0;0],[0;0;0],irf_axes(:,1),irf_axes(:,2),irf_axes(:,3),'Color',irf_color,'Linewidth',2,'ShowArrowHead','off');

if get(handles.checkbox9,'Value') == 0
    
    set(h7a,'Visible','off')
    
end

%Plot Earth Fixed Reference Frame axes
h8a                                         = quiver3([0;0;0],[0;0;0],[0;0;0],EFRF_AXES(:,1),EFRF_AXES(:,2),EFRF_AXES(:,3),'Color',efrf_color,'Linewidth',2,'ShowArrowHead','off');

if get(handles.checkbox10,'Value') == 0
    
    set(h8a,'Visible','off')
    
end

view(-24,48)
axis equal
axis off
grid off
box off
rotate3d on
hold off

%Plot to Earth Fixed Reference Frame
axes(handles.axes2)

%Plot Earth
h1b                                         = surf(X_EARTH,Y_EARTH,Z_EARTH,TOPOGRAPHY);

if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
    
    demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
    shading interp
    
elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
    
    set(h1b,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color);
    
elseif get(handles.checkbox2,'Value') == 0
    
    set(h1b,'Visible','off')
    
end

hold on

%Plot day night map
h9b                                         = surf(X_t,Y_t,Z_t,day_night_map);

shading interp
set(h9b,'Visible','off')
alpha(h9b,0.5)

%Plot coastline
h2b                                         = plot3(X_COAST,Y_COAST,Z_COAST,'Color',coastline_color);

if get(handles.checkbox4,'Value') == 0
    
    set(h2b,'Visible','off')
    
end

%Plot satelite orbit
h3b                                         = plot3(X_sat,Y_sat,Z_sat,'Color',orbit_color,'Linewidth',2);

if get(handles.checkbox6,'Value') == 0
    
    set(h3b,'Visible','off')
    
end

%Plot satellite
h4b                                         = scatter3(X_sat(1,1),Y_sat(1,1),Z_sat(1,1),100,'filled','MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);

if get(handles.checkbox5,'Value') == 0
    
    set(h4b,'Visible','off')
    
end

%Plot satellite velocity
h5b                                         = quiver3(X_sat(1,1),Y_sat(1,1),Z_sat(1,1),V_X_sat(1,1)*1000,V_Y_sat(1,1)*1000,V_Z_sat(1,1)*1000,'Color',velocity_color,'Linewidth',2);

if get(handles.checkbox7,'Value') == 0
    
    set(h5b,'Visible','off')
    
end

%Plot satellite radial distance
h6b                                         = plot3([0 X_sat(1,1)],[0 Y_sat(1,1)],[0 Z_sat(1,1)],'Color',radial_color);

if get(handles.checkbox8,'Value') == 0
    
    set(h6b,'Visible','off')
    
end

%Plot Inertial Reference Frame axes
h7b                                         = quiver3([0;0;0],[0;0;0],[0;0;0],IRF_AXES(:,1),IRF_AXES(:,2),IRF_AXES(:,3),'Color',irf_color,'Linewidth',2,'ShowArrowHead','off');

if get(handles.checkbox9,'Value') == 0
    
    set(h7b,'Visible','off')
    
end

%Plot Earth Fixed Reference Frame axes
h8b                                         = quiver3([0;0;0],[0;0;0],[0;0;0],efrf_axes(:,1),efrf_axes(:,2),efrf_axes(:,3),'Color',efrf_color,'Linewidth',2,'ShowArrowHead','off');

if get(handles.checkbox10,'Value') == 0
    
    set(h8b,'Visible','off')
    
end

view(-24,48)
axis equal
axis off
grid off
box off
rotate3d on
hold off

%Plot to Earth Fixed Reference Frame (Ground Tracks)
axes(handles.axes3)

%Plot Earth
h1c                                         = pcolor(LON,LAT,TOPOGRAPHY);

if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
    
    demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
    shading interp
    
elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
    
    set(h1c,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color);
    
elseif get(handles.checkbox2,'Value') == 0
    
    set(h1c,'Visible','off')
    
end

hold on

%Plot day night map
h9c                                         = pcolor(LONGITUDE_t,LATITUDE_t,day_night_map);

shading interp
set(h9c,'Visible','off')
alpha(h9c,0.5)

%Plot coastline
h2c                                         = plot(long_corr,lat_corr,'Color',coastline_color);

if get(handles.checkbox4,'Value') == 0
    
    set(h2c,'Visible','off')
    
end

%Plot satellite orbit
h3c                                         = plot(lamda_sat_corr,phi_sat_corr,'Color',orbit_color,'Linewidth',2);

if get(handles.checkbox6,'Value') == 0
    
    set(h3c,'Visible','off')
    
end

%Plot satellite
h4c                                         = scatter(lamda_sat(1,1),phi_sat(1,1),100,'filled','MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);
set(h4c,'ZData',1);

if get(handles.checkbox5,'Value') == 0
    
    set(h4c,'Visible','off')
    
end


axis equal
axis off
grid off
box off
rotate3d on
hold off

%Set axis auto or manual
if get(handles.checkbox1,'Value') == 1
    
    axis(handles.axes1,'equal')
    axis(handles.axes2,'equal')
    
elseif get(handles.checkbox1,'Value') == 0
    
    set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
    set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
end

addlistener(handles.slider1,'ContinuousValueChange',@(hObject, event) slider1_Callback(hObject, eventdata, handles));
addlistener(handles.slider2,'ContinuousValueChange',@(hObject, event) slider2_Callback(hObject, eventdata, handles));
addlistener(handles.slider3,'ContinuousValueChange',@(hObject, event) slider3_Callback(hObject, eventdata, handles));
addlistener(handles.slider4,'ContinuousValueChange',@(hObject, event) slider4_Callback(hObject, eventdata, handles));
addlistener(handles.slider5,'ContinuousValueChange',@(hObject, event) slider5_Callback(hObject, eventdata, handles));
addlistener(handles.slider6,'ContinuousValueChange',@(hObject, event) slider6_Callback(hObject, eventdata, handles));
addlistener(handles.slider7,'ContinuousValueChange',@(hObject, event) slider7_Callback(hObject, eventdata, handles));
addlistener(handles.slider8,'ContinuousValueChange',@(hObject, event) slider8_Callback(hObject, eventdata, handles));
addlistener(handles.slider9,'ContinuousValueChange',@(hObject, event) slider9_Callback(hObject, eventdata, handles));


% --- Executes during object deletion, before destroying properties.
function satlab_DeleteFcn(~, ~, ~)
% hObject    handle to satlab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Declare global variables
global real_time_timer current_time_timer

%Stop and delete all active timers
stop(real_time_timer)
stop(current_time_timer)

delete(real_time_timer);
delete(current_time_timer);


% --- Outputs from this function are returned to the command line.
function varargout = satlab_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

if (str2double(get(handles.edit1,'String')) >= 6371000) && (str2double(get(handles.edit1,'String')) <= 100000000)
    
    set(handles.slider1,'Value',str2double(get(handles.edit1,'String')))
    slider1_Callback(hObject, eventdata, handles)
    
else
    
    errordlg('Please select a semi-major axis value between 6371000 and 100000000 meters.')
    
end


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


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

if (str2double(get(handles.edit2,'String')) >= 0) && (str2double(get(handles.edit2,'String')) <= 1)
    
    set(handles.slider2,'Value',str2double(get(handles.edit2,'String')))
    slider2_Callback(hObject, eventdata, handles)
    
else
    
    errordlg('Eccentricity should be between 0 and 1.')
    
end


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


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

if (str2double(get(handles.edit3,'String')) >= 0) && (str2double(get(handles.edit3,'String')) <= 360)
    
    set(handles.slider3,'Value',str2double(get(handles.edit3,'String')))
    slider3_Callback(hObject, eventdata, handles)
    
else
    
    errordlg('Argument of perigee should be between 0 and 360 degrees.')
    
end


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


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

if (str2double(get(handles.edit4,'String')) >= 0) && (str2double(get(handles.edit4,'String')) <= 360)
    
    set(handles.slider4,'Value',str2double(get(handles.edit4,'String')))
    slider4_Callback(hObject, eventdata, handles)
    
else
    
    errordlg('Right ascension of ascending node should be between 0 and 360 degrees.')
    
end


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


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

if (str2double(get(handles.edit5,'String')) >= 0) && (str2double(get(handles.edit5,'String')) <= 180)
    
    set(handles.slider5,'Value',str2double(get(handles.edit5,'String')))
    slider5_Callback(hObject, eventdata, handles)
    
else
    
    errordlg('Inclination should be between 0 and 180 degrees.')
    
end


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


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

if (str2double(get(handles.edit6,'String')) >= 0) && (str2double(get(handles.edit6,'String')) <= 360)
    
    set(handles.slider6,'Value',str2double(get(handles.edit6,'String')))
    slider6_Callback(hObject, eventdata, handles)
    
else
    
    errordlg('Mean anomaly should be between 0 and 360 degrees.')
    
end


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


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

if (str2double(get(handles.edit7,'String')) >= 0) && (str2double(get(handles.edit7,'String')) <= 360)
    
    set(handles.slider7,'Value',str2double(get(handles.edit7,'String')))
    slider7_Callback(hObject, eventdata, handles)
    
else
    
    errordlg('Earth rotation angle should be between 0 and 360 degrees.')
    
end


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


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double

if (str2double(get(handles.edit8,'String')) >= 1) && (str2double(get(handles.edit8,'String')) <= 100)
    
    set(handles.slider8,'Value',str2double(get(handles.edit8,'String')))
    slider8_Callback(hObject, eventdata, handles)
    
else
    
    errordlg('Please select a value between 1 and 100 periods.')
    
end


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


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double

if (str2double(get(handles.edit9,'String')) >= 4000000) && (str2double(get(handles.edit9,'String')) <= 200000000)
    
    set(handles.slider9,'Value',str2double(get(handles.edit9,'String')))
    slider9_Callback(hObject, eventdata, handles)
    
else
    
    errordlg('Please select a value between 4000000 and 200000000 meters.')
    
end


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


function edit13_Callback(~, ~, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double

if (str2double(get(handles.edit13,'String')) > 12) || (str2double(get(handles.edit13,'String')) < -12)
        
    errordlg('Please select a value between -12 and 12 hours.')
   
else
    
    show_current_time(0, 0, handles)
    
end


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


% --- Executes on slider movement.
function slider1_Callback(~, ~, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%Declare global variables
global GM period_Earth omega_Earth O a e w W i M theta_z T_num T t era
global h3a h4a h5a h6a h3b h4b h5b h6b h3c h4c
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat phi_sat phi_sat_corr lamda_sat lamda_sat_corr

set(handles.edit1,'String',get(handles.slider1,'Value'))

%Propagate the orbit
a                                           = str2double(get(handles.edit1,'String'));
T                                           = 2*pi*sqrt(a^3/GM)/period_Earth;
t                                           = (0:T/360:T_num*T+T/360)';
era                                         = theta_z + omega_Earth*t*period_Earth;

%Inertial Reference Frame
[x_sat,y_sat,z_sat,v_x_sat,v_y_sat,v_z_sat] = f_kepl2svec(GM,a,e,w,W,i,M,0,t);

%Initialize variables
X_sat                                       = zeros(size(x_sat));
Y_sat                                       = zeros(size(y_sat));
Z_sat                                       = zeros(size(z_sat));

V_X_sat                                     = zeros(size(v_x_sat));
V_Y_sat                                     = zeros(size(v_y_sat));
V_Z_sat                                     = zeros(size(v_z_sat));

%Earth Fixed Reference Frame
for ii = 1:size(x_sat,1)
    
    pos_vec                                 = [x_sat(ii,1);y_sat(ii,1);z_sat(ii,1)];
    vel_vec                                 = [v_x_sat(ii,1);v_y_sat(ii,1);v_z_sat(ii,1)];
    
    POS_VEC                                 = f_rotz(era(ii,1))*pos_vec;
    VEL_VEC                                 = f_rotz(era(ii,1))*vel_vec - O*POS_VEC;
    
    X_sat(ii,1)                             = POS_VEC(1,1);
    Y_sat(ii,1)                             = POS_VEC(2,1);
    Z_sat(ii,1)                             = POS_VEC(3,1);
    
    V_X_sat(ii,1)                           = VEL_VEC(1,1);
    V_Y_sat(ii,1)                           = VEL_VEC(2,1);
    V_Z_sat(ii,1)                           = VEL_VEC(3,1);
    
end

%Earth Fixed Reference Frame (Map Projection)
[phi_sat,lamda_sat,~]                       = f_cart2geo(X_sat,Y_sat,Z_sat,'GRS 1980');

%Fix east-west lines bug
[phi_sat_corr,lamda_sat_corr]               = f_fix_lines(phi_sat,lamda_sat,200);

%Plot to Inertial Reference Frame
set(h3a,'XData',x_sat(1:362),'YData',y_sat(1:362),'ZData',z_sat(1:362));
set(h4a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1));
set(h5a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1),'UData',v_x_sat(1,1)*1000,'VData',v_y_sat(1,1)*1000,'WData',v_z_sat(1,1)*1000);
set(h6a,'XData',[0 x_sat(1,1)],'YData',[0 y_sat(1,1)],'ZData',[0 z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame
set(h3b,'XData',X_sat,'YData',Y_sat,'ZData',Z_sat);
set(h4b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1));
set(h5b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1),'UData',V_X_sat(1,1)*1000,'VData',V_Y_sat(1,1)*1000,'WData',V_Z_sat(1,1)*1000);
set(h6b,'XData',[0 X_sat(1,1)],'YData',[0 Y_sat(1,1)],'ZData',[0 Z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame (Map Projection)
set(h3c,'XData',lamda_sat_corr,'YData',phi_sat_corr);
set(h4c,'XData',lamda_sat(1,1),'YData',phi_sat(1,1));

%Set axis auto or manual
if get(handles.checkbox1,'Value') == 1
    
    axis(handles.axes1,'equal')
    axis(handles.axes2,'equal')
    
elseif get(handles.checkbox1,'Value') == 0
    
    set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
    set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
end


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, ~, ~)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(~, ~, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%Declare global variables
global GM O a e w W i M t era
global h3a h4a h5a h6a h3b h4b h5b h6b h3c h4c
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat  phi_sat phi_sat_corr lamda_sat lamda_sat_corr

set(handles.edit2,'String',get(handles.slider2,'Value'))

%Propagate the orbit
e                                           = str2double(get(handles.edit2,'String'));

%Inertial Reference Frame
[x_sat,y_sat,z_sat,v_x_sat,v_y_sat,v_z_sat] = f_kepl2svec(GM,a,e,w,W,i,M,0,t);

%Initialize variables
X_sat                                       = zeros(size(x_sat));
Y_sat                                       = zeros(size(y_sat));
Z_sat                                       = zeros(size(z_sat));

V_X_sat                                     = zeros(size(v_x_sat));
V_Y_sat                                     = zeros(size(v_y_sat));
V_Z_sat                                     = zeros(size(v_z_sat));

%Earth Fixed Reference Frame
for ii = 1:size(x_sat,1)
    
    pos_vec                                 = [x_sat(ii,1);y_sat(ii,1);z_sat(ii,1)];
    vel_vec                                 = [v_x_sat(ii,1);v_y_sat(ii,1);v_z_sat(ii,1)];
    
    POS_VEC                                 = f_rotz(era(ii,1))*pos_vec;
    VEL_VEC                                 = f_rotz(era(ii,1))*vel_vec - O*POS_VEC;
    
    X_sat(ii,1)                             = POS_VEC(1,1);
    Y_sat(ii,1)                             = POS_VEC(2,1);
    Z_sat(ii,1)                             = POS_VEC(3,1);
    
    V_X_sat(ii,1)                           = VEL_VEC(1,1);
    V_Y_sat(ii,1)                           = VEL_VEC(2,1);
    V_Z_sat(ii,1)                           = VEL_VEC(3,1);
    
end

%Earth Fixed Reference Frame (Map Projection)
[phi_sat,lamda_sat,~]                       = f_cart2geo(X_sat,Y_sat,Z_sat,'GRS 1980');

%Fix east-west lines bug
[phi_sat_corr,lamda_sat_corr]               = f_fix_lines(phi_sat,lamda_sat,200);

%Plot to Inertial Reference Frame
set(h3a,'XData',x_sat(1:362),'YData',y_sat(1:362),'ZData',z_sat(1:362));
set(h4a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1));
set(h5a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1),'UData',v_x_sat(1,1)*1000,'VData',v_y_sat(1,1)*1000,'WData',v_z_sat(1,1)*1000);
set(h6a,'XData',[0 x_sat(1,1)],'YData',[0 y_sat(1,1)],'ZData',[0 z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame
set(h3b,'XData',X_sat,'YData',Y_sat,'ZData',Z_sat);
set(h4b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1));
set(h5b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1),'UData',V_X_sat(1,1)*1000,'VData',V_Y_sat(1,1)*1000,'WData',V_Z_sat(1,1)*1000);
set(h6b,'XData',[0 X_sat(1,1)],'YData',[0 Y_sat(1,1)],'ZData',[0 Z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame (Map Projection)
set(h3c,'XData',lamda_sat_corr,'YData',phi_sat_corr);
set(h4c,'XData',lamda_sat(1,1),'YData',phi_sat(1,1));

%Set axis auto or manual
if get(handles.checkbox1,'Value') == 1
    
    axis(handles.axes1,'equal')
    axis(handles.axes2,'equal')
    
elseif get(handles.checkbox1,'Value') == 0
    
    set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
    set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
end


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, ~, ~)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(~, ~, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%Declare global variables
global GM O a e w W i M t era
global h3a h4a h5a h6a h3b h4b h5b h6b h3c h4c
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat phi_sat phi_sat_corr lamda_sat lamda_sat_corr

set(handles.edit3,'String',get(handles.slider3,'Value'))

%Propagate the orbit
w                                           = str2double(get(handles.edit3,'String'));

%Inertial Reference Frame
[x_sat,y_sat,z_sat,v_x_sat,v_y_sat,v_z_sat] = f_kepl2svec(GM,a,e,w,W,i,M,0,t);

%Initialize variables
X_sat                                       = zeros(size(x_sat));
Y_sat                                       = zeros(size(y_sat));
Z_sat                                       = zeros(size(z_sat));

V_X_sat                                     = zeros(size(v_x_sat));
V_Y_sat                                     = zeros(size(v_y_sat));
V_Z_sat                                     = zeros(size(v_z_sat));

%Earth Fixed Reference Frame
for ii = 1:size(x_sat,1)
    
    pos_vec                                 = [x_sat(ii,1);y_sat(ii,1);z_sat(ii,1)];
    vel_vec                                 = [v_x_sat(ii,1);v_y_sat(ii,1);v_z_sat(ii,1)];
    
    POS_VEC                                 = f_rotz(era(ii,1))*pos_vec;
    VEL_VEC                                 = f_rotz(era(ii,1))*vel_vec - O*POS_VEC;
    
    X_sat(ii,1)                             = POS_VEC(1,1);
    Y_sat(ii,1)                             = POS_VEC(2,1);
    Z_sat(ii,1)                             = POS_VEC(3,1);
    
    V_X_sat(ii,1)                           = VEL_VEC(1,1);
    V_Y_sat(ii,1)                           = VEL_VEC(2,1);
    V_Z_sat(ii,1)                           = VEL_VEC(3,1);
    
end

%Earth Fixed Reference Frame (Map Projection)
[phi_sat,lamda_sat,~]                       = f_cart2geo(X_sat,Y_sat,Z_sat,'GRS 1980');

%Fix east-west lines bug
[phi_sat_corr,lamda_sat_corr]               = f_fix_lines(phi_sat,lamda_sat,200);

%Plot to Inertial Reference Frame
set(h3a,'XData',x_sat(1:362),'YData',y_sat(1:362),'ZData',z_sat(1:362));
set(h4a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1));
set(h5a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1),'UData',v_x_sat(1,1)*1000,'VData',v_y_sat(1,1)*1000,'WData',v_z_sat(1,1)*1000);
set(h6a,'XData',[0 x_sat(1,1)],'YData',[0 y_sat(1,1)],'ZData',[0 z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame
set(h3b,'XData',X_sat,'YData',Y_sat,'ZData',Z_sat);
set(h4b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1));
set(h5b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1),'UData',V_X_sat(1,1)*1000,'VData',V_Y_sat(1,1)*1000,'WData',V_Z_sat(1,1)*1000);
set(h6b,'XData',[0 X_sat(1,1)],'YData',[0 Y_sat(1,1)],'ZData',[0 Z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame (Map Projection)
set(h3c,'XData',lamda_sat_corr,'YData',phi_sat_corr);
set(h4c,'XData',lamda_sat(1,1),'YData',phi_sat(1,1));

%Set axis auto or manual
if get(handles.checkbox1,'Value') == 1
    
    axis(handles.axes1,'equal')
    axis(handles.axes2,'equal')
    
elseif get(handles.checkbox1,'Value') == 0
    
    set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
    set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
end


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, ~, ~)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(~, ~, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%Declare global variables
global GM O a e w W i M t era
global h3a h4a h5a h6a h3b h4b h5b h6b h3c h4c
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat phi_sat phi_sat_corr lamda_sat lamda_sat_corr

set(handles.edit4,'String',get(handles.slider4,'Value'))

%Propagate the orbit
W                                           = str2double(get(handles.edit4,'String'));

%Inertial Reference Frame
[x_sat,y_sat,z_sat,v_x_sat,v_y_sat,v_z_sat] = f_kepl2svec(GM,a,e,w,W,i,M,0,t);

%Initialize variables
X_sat                                       = zeros(size(x_sat));
Y_sat                                       = zeros(size(y_sat));
Z_sat                                       = zeros(size(z_sat));

V_X_sat                                     = zeros(size(v_x_sat));
V_Y_sat                                     = zeros(size(v_y_sat));
V_Z_sat                                     = zeros(size(v_z_sat));

%Earth Fixed Reference Frame
for ii = 1:size(x_sat,1)
    
    pos_vec                                 = [x_sat(ii,1);y_sat(ii,1);z_sat(ii,1)];
    vel_vec                                 = [v_x_sat(ii,1);v_y_sat(ii,1);v_z_sat(ii,1)];
    
    POS_VEC                                 = f_rotz(era(ii,1))*pos_vec;
    VEL_VEC                                 = f_rotz(era(ii,1))*vel_vec - O*POS_VEC;
    
    X_sat(ii,1)                             = POS_VEC(1,1);
    Y_sat(ii,1)                             = POS_VEC(2,1);
    Z_sat(ii,1)                             = POS_VEC(3,1);
    
    V_X_sat(ii,1)                           = VEL_VEC(1,1);
    V_Y_sat(ii,1)                           = VEL_VEC(2,1);
    V_Z_sat(ii,1)                           = VEL_VEC(3,1);
    
end

%Earth Fixed Reference Frame (Map Projection)
[phi_sat,lamda_sat,~]                       = f_cart2geo(X_sat,Y_sat,Z_sat,'GRS 1980');

%Fix east-west lines bug
[phi_sat_corr,lamda_sat_corr]               = f_fix_lines(phi_sat,lamda_sat,200);

%Plot to Inertial Reference Frame
set(h3a,'XData',x_sat(1:362),'YData',y_sat(1:362),'ZData',z_sat(1:362));
set(h4a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1));
set(h5a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1),'UData',v_x_sat(1,1)*1000,'VData',v_y_sat(1,1)*1000,'WData',v_z_sat(1,1)*1000);
set(h6a,'XData',[0 x_sat(1,1)],'YData',[0 y_sat(1,1)],'ZData',[0 z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame
set(h3b,'XData',X_sat,'YData',Y_sat,'ZData',Z_sat);
set(h4b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1));
set(h5b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1),'UData',V_X_sat(1,1)*1000,'VData',V_Y_sat(1,1)*1000,'WData',V_Z_sat(1,1)*1000);
set(h6b,'XData',[0 X_sat(1,1)],'YData',[0 Y_sat(1,1)],'ZData',[0 Z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame (Map Projection)
set(h3c,'XData',lamda_sat_corr,'YData',phi_sat_corr);
set(h4c,'XData',lamda_sat(1,1),'YData',phi_sat(1,1));

%Set axis auto or manual
if get(handles.checkbox1,'Value') == 1
    
    axis(handles.axes1,'equal')
    axis(handles.axes2,'equal')
    
elseif get(handles.checkbox1,'Value') == 0
    
    set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
    set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
end


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, ~, ~)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(~, ~, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%Declare global variables
global GM O a e w W i M t era
global h3a h4a h5a h6a h3b h4b h5b h6b h3c h4c
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat phi_sat phi_sat_corr lamda_sat lamda_sat_corr

set(handles.edit5,'String',get(handles.slider5,'Value'))

%Propagate the orbit
i                                           = str2double(get(handles.edit5,'String'));

%Inertial Reference Frame
[x_sat,y_sat,z_sat,v_x_sat,v_y_sat,v_z_sat] = f_kepl2svec(GM,a,e,w,W,i,M,0,t);

%Initialize variables
X_sat                                       = zeros(size(x_sat));
Y_sat                                       = zeros(size(y_sat));
Z_sat                                       = zeros(size(z_sat));

V_X_sat                                     = zeros(size(v_x_sat));
V_Y_sat                                     = zeros(size(v_y_sat));
V_Z_sat                                     = zeros(size(v_z_sat));

%Earth Fixed Reference Frame
for ii = 1:size(x_sat,1)
    
    pos_vec                                 = [x_sat(ii,1);y_sat(ii,1);z_sat(ii,1)];
    vel_vec                                 = [v_x_sat(ii,1);v_y_sat(ii,1);v_z_sat(ii,1)];
    
    POS_VEC                                 = f_rotz(era(ii,1))*pos_vec;
    VEL_VEC                                 = f_rotz(era(ii,1))*vel_vec - O*POS_VEC;
    
    X_sat(ii,1)                             = POS_VEC(1,1);
    Y_sat(ii,1)                             = POS_VEC(2,1);
    Z_sat(ii,1)                             = POS_VEC(3,1);
    
    V_X_sat(ii,1)                           = VEL_VEC(1,1);
    V_Y_sat(ii,1)                           = VEL_VEC(2,1);
    V_Z_sat(ii,1)                           = VEL_VEC(3,1);
    
end

%Earth Fixed Reference Frame (Map Projection)
[phi_sat,lamda_sat,~]                       = f_cart2geo(X_sat,Y_sat,Z_sat,'GRS 1980');

%Fix east-west lines bug
[phi_sat_corr,lamda_sat_corr]               = f_fix_lines(phi_sat,lamda_sat,200);

%Plot to Inertial Reference Frame
set(h3a,'XData',x_sat(1:362),'YData',y_sat(1:362),'ZData',z_sat(1:362));
set(h4a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1));
set(h5a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1),'UData',v_x_sat(1,1)*1000,'VData',v_y_sat(1,1)*1000,'WData',v_z_sat(1,1)*1000);
set(h6a,'XData',[0 x_sat(1,1)],'YData',[0 y_sat(1,1)],'ZData',[0 z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame
set(h3b,'XData',X_sat,'YData',Y_sat,'ZData',Z_sat);
set(h4b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1));
set(h5b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1),'UData',V_X_sat(1,1)*1000,'VData',V_Y_sat(1,1)*1000,'WData',V_Z_sat(1,1)*1000);
set(h6b,'XData',[0 X_sat(1,1)],'YData',[0 Y_sat(1,1)],'ZData',[0 Z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame (Map Projection)
set(h3c,'XData',lamda_sat_corr,'YData',phi_sat_corr);
set(h4c,'XData',lamda_sat(1,1),'YData',phi_sat(1,1));

%Set axis auto or manual
if get(handles.checkbox1,'Value') == 1
    
    axis(handles.axes1,'equal')
    axis(handles.axes2,'equal')
    
elseif get(handles.checkbox1,'Value') == 0
    
    set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
    set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
end


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, ~, ~)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(~, ~, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%Declare global variables
global GM O a e w W i M t era
global h3a h4a h5a h6a h3b h4b h5b h6b h3c h4c
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat phi_sat phi_sat_corr lamda_sat lamda_sat_corr

set(handles.edit6,'String',get(handles.slider6,'Value'))

%Propagate the orbit
M                                           = str2double(get(handles.edit6,'String'));

%Inertial Reference Frame
[x_sat,y_sat,z_sat,v_x_sat,v_y_sat,v_z_sat] = f_kepl2svec(GM,a,e,w,W,i,M,0,t);

%Initialize variables
X_sat                                       = zeros(size(x_sat));
Y_sat                                       = zeros(size(y_sat));
Z_sat                                       = zeros(size(z_sat));

V_X_sat                                     = zeros(size(v_x_sat));
V_Y_sat                                     = zeros(size(v_y_sat));
V_Z_sat                                     = zeros(size(v_z_sat));

%Earth Fixed Reference Frame
for ii = 1:size(x_sat,1)
    
    pos_vec                                 = [x_sat(ii,1);y_sat(ii,1);z_sat(ii,1)];
    vel_vec                                 = [v_x_sat(ii,1);v_y_sat(ii,1);v_z_sat(ii,1)];
    
    POS_VEC                                 = f_rotz(era(ii,1))*pos_vec;
    VEL_VEC                                 = f_rotz(era(ii,1))*vel_vec - O*POS_VEC;
    
    X_sat(ii,1)                             = POS_VEC(1,1);
    Y_sat(ii,1)                             = POS_VEC(2,1);
    Z_sat(ii,1)                             = POS_VEC(3,1);
    
    V_X_sat(ii,1)                           = VEL_VEC(1,1);
    V_Y_sat(ii,1)                           = VEL_VEC(2,1);
    V_Z_sat(ii,1)                           = VEL_VEC(3,1);
    
end

%Earth Fixed Reference Frame (Map Projection)
[phi_sat,lamda_sat,~]                       = f_cart2geo(X_sat,Y_sat,Z_sat,'GRS 1980');

%Fix east-west lines bug
[phi_sat_corr,lamda_sat_corr]               = f_fix_lines(phi_sat,lamda_sat,200);

%Plot to Inertial Reference Frame
set(h3a,'XData',x_sat(1:362),'YData',y_sat(1:362),'ZData',z_sat(1:362));
set(h4a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1));
set(h5a,'XData',x_sat(1,1),'YData',y_sat(1,1),'ZData',z_sat(1,1),'UData',v_x_sat(1,1)*1000,'VData',v_y_sat(1,1)*1000,'WData',v_z_sat(1,1)*1000);
set(h6a,'XData',[0 x_sat(1,1)],'YData',[0 y_sat(1,1)],'ZData',[0 z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame
set(h3b,'XData',X_sat,'YData',Y_sat,'ZData',Z_sat);
set(h4b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1));
set(h5b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1),'UData',V_X_sat(1,1)*1000,'VData',V_Y_sat(1,1)*1000,'WData',V_Z_sat(1,1)*1000);
set(h6b,'XData',[0 X_sat(1,1)],'YData',[0 Y_sat(1,1)],'ZData',[0 Z_sat(1,1)]);

%Plot to Earth Fixed Reference Frame (Map Projection)
set(h3c,'XData',lamda_sat_corr,'YData',phi_sat_corr);
set(h4c,'XData',lamda_sat(1,1),'YData',phi_sat(1,1));

%Set axis auto or manual
if get(handles.checkbox1,'Value') == 1
    
    axis(handles.axes1,'equal')
    axis(handles.axes2,'equal')
    
elseif get(handles.checkbox1,'Value') == 0
    
    set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
    set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
end


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, ~, ~)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider7_Callback(~, ~, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%Declare global variables
global POS_EARTH POS_COAST row col irf_axes efrf_axes
global GM period_Earth omega_Earth O a e w W i M theta_z t era
global h1a h2a h8a h3b h4b h5b h6b h7b h3c h4c
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat phi_sat phi_sat_corr lamda_sat lamda_sat_corr

set(handles.edit7,'String',get(handles.slider7,'Value'))

%Inertial Reference Frame
[x_sat,y_sat,z_sat,v_x_sat,v_y_sat,v_z_sat] = f_kepl2svec(GM,a,e,w,W,i,M,0,t);

theta_z                                     = get(handles.slider7,'Value')*pi/180;
era                                         = theta_z + omega_Earth*t*period_Earth;

%Initialize variables
X_sat                                       = zeros(size(x_sat));
Y_sat                                       = zeros(size(y_sat));
Z_sat                                       = zeros(size(z_sat));

V_X_sat                                     = zeros(size(v_x_sat));
V_Y_sat                                     = zeros(size(v_y_sat));
V_Z_sat                                     = zeros(size(v_z_sat));

%Earth Fixed Reference Frame
for ii = 1:size(x_sat,1)
    
    pos_vec                                 = [x_sat(ii,1);y_sat(ii,1);z_sat(ii,1)];
    vel_vec                                 = [v_x_sat(ii,1);v_y_sat(ii,1);v_z_sat(ii,1)];
    
    POS_VEC                                 = f_rotz(era(ii,1))*pos_vec;
    VEL_VEC                                 = f_rotz(era(ii,1))*vel_vec - O*POS_VEC;
    
    X_sat(ii,1)                             = POS_VEC(1,1);
    Y_sat(ii,1)                             = POS_VEC(2,1);
    Z_sat(ii,1)                             = POS_VEC(3,1);
    
    V_X_sat(ii,1)                           = VEL_VEC(1,1);
    V_Y_sat(ii,1)                           = VEL_VEC(2,1);
    V_Z_sat(ii,1)                           = VEL_VEC(3,1);
    
end

POS_EARTH_NEW                               = f_rotz(-theta_z)*POS_EARTH;
POS_COAST_NEW                               = f_rotz(-theta_z)*POS_COAST;

X_EARTH_NEW                                 = reshape(POS_EARTH_NEW(1,:),row,col);
Y_EARTH_NEW                                 = reshape(POS_EARTH_NEW(2,:),row,col);
Z_EARTH_NEW                                 = reshape(POS_EARTH_NEW(3,:),row,col);

X_COAST_NEW                                 = POS_COAST_NEW(1,:)';
Y_COAST_NEW                                 = POS_COAST_NEW(2,:)';
Z_COAST_NEW                                 = POS_COAST_NEW(3,:)';

IRF_AXES                                    = f_rotz(-theta_z)*irf_axes;
EFRF_AXES                                   = f_rotz(theta_z)*efrf_axes;

%Earth Fixed Reference Frame (Map Projection)
[phi_sat,lamda_sat,~]                       = f_cart2geo(X_sat,Y_sat,Z_sat,'GRS 1980');

%Fix east-west lines bug
[phi_sat_corr,lamda_sat_corr]               = f_fix_lines(phi_sat,lamda_sat,200);

%Plot to Inertial Reference Frame
set(h1a,'XData',X_EARTH_NEW,'YData',Y_EARTH_NEW,'ZData',Z_EARTH_NEW);
set(h2a,'XData',X_COAST_NEW,'YData',Y_COAST_NEW,'ZData',Z_COAST_NEW);
set(h8a,'XData',[0;0;0],'YData',[0;0;0],'ZData',[0;0;0],'UData',EFRF_AXES(:,1),'VData',EFRF_AXES(:,2),'WData',EFRF_AXES(:,3));

%Plot to Earth Fixed Reference Frame
set(h3b,'XData',X_sat,'YData',Y_sat,'ZData',Z_sat);
set(h4b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1));
set(h5b,'XData',X_sat(1,1),'YData',Y_sat(1,1),'ZData',Z_sat(1,1),'UData',V_X_sat(1,1)*1000,'VData',V_Y_sat(1,1)*1000,'WData',V_Z_sat(1,1)*1000);
set(h6b,'XData',[0 X_sat(1,1)],'YData',[0 Y_sat(1,1)],'ZData',[0 Z_sat(1,1)]);
set(h7b,'XData',[0;0;0],'YData',[0;0;0],'ZData',[0;0;0],'UData',IRF_AXES(:,1),'VData',IRF_AXES(:,2),'WData',IRF_AXES(:,3));

%Plot to Earth Fixed Reference Frame (Map Projection)
set(h3c,'XData',lamda_sat_corr,'YData',phi_sat_corr);
set(h4c,'XData',lamda_sat(1,1),'YData',phi_sat(1,1));

%Set axis auto or manual
if get(handles.checkbox1,'Value') == 1
    
    axis(handles.axes1,'equal')
    axis(handles.axes2,'equal')
    
elseif get(handles.checkbox1,'Value') == 0
    
    set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
    set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
end


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, ~, ~)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider8_Callback(~, ~, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%Declare global variables
global GM period_Earth omega_Earth O a e w W i M theta_z T_num T t era
global h3a h3b h3c
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat phi_sat phi_sat_corr lamda_sat lamda_sat_corr

set(handles.edit8,'String',get(handles.slider8,'Value'))

T_num                                       = str2double(get(handles.edit8,'String'));
t                                           = (0:T/360:T_num*T+T/360)';
era                                         = theta_z + omega_Earth*t*period_Earth;

%Inertial Reference Frame
[x_sat,y_sat,z_sat,v_x_sat,v_y_sat,v_z_sat] = f_kepl2svec(GM,a,e,w,W,i,M,0,t);

%Initialize variables
X_sat                                       = zeros(size(x_sat));
Y_sat                                       = zeros(size(y_sat));
Z_sat                                       = zeros(size(z_sat));

V_X_sat                                     = zeros(size(v_x_sat));
V_Y_sat                                     = zeros(size(v_y_sat));
V_Z_sat                                     = zeros(size(v_z_sat));

%Earth Fixed Reference Frame
for ii = 1:size(x_sat,1)
    
    pos_vec                                 = [x_sat(ii,1);y_sat(ii,1);z_sat(ii,1)];
    vel_vec                                 = [v_x_sat(ii,1);v_y_sat(ii,1);v_z_sat(ii,1)];
    
    POS_VEC                                 = f_rotz(era(ii,1))*pos_vec;
    VEL_VEC                                 = f_rotz(era(ii,1))*vel_vec - O*POS_VEC;
    
    X_sat(ii,1)                             = POS_VEC(1,1);
    Y_sat(ii,1)                             = POS_VEC(2,1);
    Z_sat(ii,1)                             = POS_VEC(3,1);
    
    V_X_sat(ii,1)                           = VEL_VEC(1,1);
    V_Y_sat(ii,1)                           = VEL_VEC(2,1);
    V_Z_sat(ii,1)                           = VEL_VEC(3,1);
    
end

%Earth Fixed Reference Frame (Map Projection)
[phi_sat,lamda_sat,~]                       = f_cart2geo(X_sat,Y_sat,Z_sat,'GRS 1980');

%Fix east-west lines bug
[phi_sat_corr,lamda_sat_corr]               = f_fix_lines(phi_sat,lamda_sat,200);

%Plot to Inertial Reference Frame
set(h3a,'XData',x_sat(1:362),'YData',y_sat(1:362),'ZData',z_sat(1:362));

%Plot to Earth Fixed Reference Frame
set(h3b,'XData',X_sat,'YData',Y_sat,'ZData',Z_sat);

%Plot to Earth Fixed Reference Frame (Map Projection)
set(h3c,'XData',lamda_sat_corr,'YData',phi_sat_corr);

%Set axis auto or manual
if get(handles.checkbox1,'Value') == 1
    
    axis(handles.axes1,'equal')
    axis(handles.axes2,'equal')
    
elseif get(handles.checkbox1,'Value') == 0
    
    set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
    set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
end


% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, ~, ~)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider9_Callback(~, ~, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.edit9,'String',get(handles.slider9,'Value'))

set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])

set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])


% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, ~, ~)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(~, ~, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

if get(handles.checkbox1,'Value') == 0
    
    set(handles.edit9,'Enable','On')
    set(handles.slider9,'Enable','On')
    
elseif get(handles.checkbox1,'Value') == 1
    
    set(handles.edit9,'Enable','Off')
    set(handles.slider9,'Enable','Off')
    
end

%Set axis auto or manual
if get(handles.checkbox1,'Value') == 1
    
    axis(handles.axes1,'equal')
    axis(handles.axes2,'equal')
    
elseif get(handles.checkbox1,'Value') == 0
    
    set(handles.axes1,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes1,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
    set(handles.axes2,'XLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'YLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    set(handles.axes2,'ZLim',[-get(handles.slider9,'Value') get(handles.slider9,'Value')])
    
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(~, ~, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

%Declare global variables
global h1a h1b h1c TOPOGRAPHY earth_color

if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
    
    set(handles.checkbox3,'Enable','on')
    
    set(h1a,'CData',TOPOGRAPHY,'EdgeColor','none','FaceColor','interp')
    set(h1b,'CData',TOPOGRAPHY,'EdgeColor','none','FaceColor','interp')
    set(h1c,'CData',TOPOGRAPHY,'EdgeColor','none','FaceColor','interp')
    
    set(h1a,'Visible','on')
    set(h1b,'Visible','on')
    set(h1c,'Visible','on')
    
    demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
    
elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
    
    set(handles.checkbox3,'Enable','on')
    
    set(h1a,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color)
    set(h1b,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color)
    set(h1c,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color)
    
    set(h1a,'Visible','on')
    set(h1b,'Visible','on')
    set(h1c,'Visible','on')
    
elseif get(handles.checkbox2,'Value') == 0
    
    set(handles.checkbox3,'Enable','off')
    
    set(h1a,'Visible','off')
    set(h1b,'Visible','off')
    set(h1c,'Visible','off')
    
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(~, ~, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3

%Declare global variables
global h1a h1b h1c TOPOGRAPHY earth_color

if get(handles.checkbox3,'Value') == 1
    
    set(h1a,'CData',TOPOGRAPHY,'EdgeColor','none','FaceColor','interp')
    set(h1b,'CData',TOPOGRAPHY,'EdgeColor','none','FaceColor','interp')
    set(h1c,'CData',TOPOGRAPHY,'EdgeColor','none','FaceColor','interp')
    
    demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
    
elseif get(handles.checkbox3,'Value') == 0
    
    set(h1a,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color)
    set(h1b,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color)
    set(h1c,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color)
    
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(~, ~, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4

%Declare global variables
global h2a h2b h2c

if get(handles.checkbox4,'Value') == 0
    
    set(h2a,'Visible','off')
    set(h2b,'Visible','off')
    set(h2c,'Visible','off')
    
elseif get(handles.checkbox4,'Value') == 1
    
    set(h2a,'Visible','on')
    set(h2b,'Visible','on')
    set(h2c,'Visible','on')
    
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(~, ~, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5

%Declare global variables
global h4a h4b h4c

if get(handles.checkbox5,'Value') == 0
    
    set(h4a,'Visible','off')
    set(h4b,'Visible','off')
    set(h4c,'Visible','off')
    
elseif get(handles.checkbox5,'Value') == 1
    
    set(h4a,'Visible','on')
    set(h4b,'Visible','on')
    set(h4c,'Visible','on')
    
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(~, ~, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6

%Declare global variables
global h3a h3b h3c

if get(handles.checkbox6,'Value') == 0
    
    set(h3a,'Visible','off')
    set(h3b,'Visible','off')
    set(h3c,'Visible','off')
    
elseif get(handles.checkbox6,'Value') == 1
    
    set(h3a,'Visible','on')
    set(h3b,'Visible','on')
    set(h3c,'Visible','on')
    
end


% --- Executes on button press in checkbox7.
function checkbox7_Callback(~, ~, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7

%Declare global variables
global h5a h5b

if get(handles.checkbox7,'Value') == 0
    
    set(h5a,'Visible','off')
    set(h5b,'Visible','off')
    
elseif get(handles.checkbox7,'Value') == 1
    
    set(h5a,'Visible','on')
    set(h5b,'Visible','on')
    
end


% --- Executes on button press in checkbox8.
function checkbox8_Callback(~, ~, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8

%Declare global variables
global h6a h6b

if get(handles.checkbox8,'Value') == 0
    
    set(h6a,'Visible','off')
    set(h6b,'Visible','off')
    
elseif get(handles.checkbox8,'Value') == 1
    
    set(h6a,'Visible','on')
    set(h6b,'Visible','on')
    
end


% --- Executes on button press in checkbox9.
function checkbox9_Callback(~, ~, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9

%Declare global variables
global h7a h7b

if get(handles.checkbox9,'Value') == 0
    
    set(h7a,'Visible','off')
    set(h7b,'Visible','off')
    
elseif get(handles.checkbox9,'Value') == 1
    
    set(h7a,'Visible','on')
    set(h7b,'Visible','on')
    
end


% --- Executes on button press in checkbox10.
function checkbox10_Callback(~, ~, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10

%Declare global variables
global h8a h8b

if get(handles.checkbox10,'Value') == 0
    
    set(h8a,'Visible','off')
    set(h8b,'Visible','off')
    
elseif get(handles.checkbox10,'Value') == 1
    
    set(h8a,'Visible','on')
    set(h8b,'Visible','on')
    
end


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11

%Declare global variables
global GM a e w W i M theta_z
global sat_data_new sat_num real_time_timer reference_epoch current_epoch
global h9a h9b h9c

sat_num                                         = get(handles.popupmenu1,'Value');
sat_num                                         = sat_num - 1;

if get(handles.checkbox11,'Value') == 1 && sat_num ~=0
    
    set(handles.checkbox12,'Enable','on')
    
    set(handles.edit1,'Enable','off')
    set(handles.edit2,'Enable','off')
    set(handles.edit3,'Enable','off')
    set(handles.edit4,'Enable','off')
    set(handles.edit5,'Enable','off')
    set(handles.edit6,'Enable','off')
    set(handles.edit7,'Enable','off')
    
    set(handles.slider1,'Enable','off')
    set(handles.slider2,'Enable','off')
    set(handles.slider3,'Enable','off')
    set(handles.slider4,'Enable','off')
    set(handles.slider5,'Enable','off')
    set(handles.slider6,'Enable','off')
    set(handles.slider7,'Enable','off')
    
    a                                           = (GM/(sat_data_new(sat_num).mean_motion^2))^(1/3);
    e                                           = sat_data_new(sat_num).eccentricity;
    w                                           = sat_data_new(sat_num).argument_of_perigee   ;
    W                                           = sat_data_new(sat_num).right_ascension_AN;
    i                                           = sat_data_new(sat_num).inclination;
    M                                           = sat_data_new(sat_num).mean_anomaly;
    
    theta_z                                     = f_ut12lmst(0,f_gd2jd(sat_data_new(sat_num).epoch_day,1,sat_data_new(sat_num).epoch_year,0,0,0,'JD'),'JD');
    reference_epoch                             = f_gd2jd(sat_data_new(sat_num).epoch_day,1,sat_data_new(sat_num).epoch_year,0,0,0,'JD');
    
    set(handles.edit1,'String',a)
    set(handles.edit2,'String',e)
    set(handles.edit3,'String',w)
    set(handles.edit4,'String',W)
    set(handles.edit5,'String',i)
    set(handles.edit6,'String',M)
    set(handles.edit7,'String',theta_z)
    set(handles.edit10,'String',num2str(reference_epoch))
    
    if (current_epoch - reference_epoch > 2) & (get(handles.checkbox12,'Value') == 1)
        
        set(handles.text14,'String','CAUTION! Outdated orbital data.')
        set(handles.text14,'ForeGroundColor','red')
        
    elseif (current_epoch - reference_epoch <= 2) & (get(handles.checkbox12,'Value') == 1)
        
        set(handles.text14,'String','Orbital data up-to-date!')
        set(handles.text14,'ForeGroundColor','green')
        
    end
    
    edit1_Callback(hObject, eventdata, handles)
    edit2_Callback(hObject, eventdata, handles)
    edit3_Callback(hObject, eventdata, handles)
    edit4_Callback(hObject, eventdata, handles)
    edit5_Callback(hObject, eventdata, handles)
    edit6_Callback(hObject, eventdata, handles)
    edit7_Callback(hObject, eventdata, handles)
    edit8_Callback(hObject, eventdata, handles)
    
else
    
    set(handles.checkbox12,'Value',0)
    set(handles.checkbox12,'Enable','off')
    stop(real_time_timer)
    
    set(h9a,'Visible','off')
    set(h9b,'Visible','off')
    set(h9c,'Visible','off')
    
    set(handles.edit1,'Enable','on')
    set(handles.edit2,'Enable','on')
    set(handles.edit3,'Enable','on')
    set(handles.edit4,'Enable','on')
    set(handles.edit5,'Enable','on')
    set(handles.edit6,'Enable','on')
    set(handles.edit7,'Enable','on')
    set(handles.edit8,'Enable','on')
    
    set(handles.slider1,'Enable','on')
    set(handles.slider2,'Enable','on')
    set(handles.slider3,'Enable','on')
    set(handles.slider4,'Enable','on')
    set(handles.slider5,'Enable','on')
    set(handles.slider6,'Enable','on')
    set(handles.slider7,'Enable','on')
    set(handles.slider8,'Enable','on')
    
    set(handles.edit10,'String','')
    set(handles.text14,'String','')
    set(handles.text14,'ForeGroundColor','k')
    
end


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12

%Declare global variables
global real_time_timer reference_epoch current_epoch
global h9a h9b h9c

if get(handles.checkbox12,'Value') == 1 && strcmp(get(handles.checkbox12,'Enable'),'on') == 1
    
    start(real_time_timer)
    
    set(handles.edit8,'Enable','off')
    set(handles.slider8,'Enable','off')
    
    if current_epoch - reference_epoch > 2
        
        set(handles.text14,'String','CAUTION! Outdated orbital data.')
        set(handles.text14,'ForeGroundColor','red')
        
    else
        
        set(handles.text14,'String','Orbital data up-to-date!')
        set(handles.text14,'ForeGroundColor','green')
        
    end
    
else
    
    stop(real_time_timer)
    
    set(h9a,'Visible','off')
    set(h9b,'Visible','off')
    set(h9c,'Visible','off')
    
    set(handles.edit8,'Enable','on')
    set(handles.slider8,'Enable','on')
    
    set(handles.text14,'String','')
    set(handles.text14,'ForeGroundColor','k')
    
    checkbox11_Callback(hObject, eventdata, handles)
    
end


% --- Executes on button press in checkbox15.
function checkbox15_Callback(~, ~, ~)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15


% --- Executes on button press in checkbox16.
function checkbox16_Callback(~, ~, ~)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16


% --- Executes on button press in checkbox17.
function checkbox17_Callback(~, ~, ~)
% hObject    handle to checkbox17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Declare global variables
global h1a h1b h1c earth_color

earth_color     = uisetcolor;

if length(earth_color) ~= 1
    
    set(handles.pushbutton1,'BackGroundColor',earth_color)
    
    if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
        
        set(h1a,'FaceColor',earth_color);
        set(h1b,'FaceColor',earth_color);
        set(h1c,'FaceColor',earth_color);
        
    end
    
elseif length(earth_color) == 1
    
    earth_color = get(handles.pushbutton1,'BackGroundColor');
    
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(~, ~, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Declare global variables
global h2a h2b h2c coastline_color

coastline_color     = uisetcolor;

if length(coastline_color) ~= 1
    
    set(handles.pushbutton2,'BackGroundColor',coastline_color)
    
    set(h2a,'Color',coastline_color);
    set(h2b,'Color',coastline_color);
    set(h2c,'Color',coastline_color);
    
elseif length(coastline_color) == 1
    
    coastline_color = get(handles.pushbutton2,'BackGroundColor');
    
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(~, ~, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Declare global variables
global h4a h4b h4c satellite_color

satellite_color     = uisetcolor;

if length(satellite_color) ~= 1
    
    set(handles.pushbutton3,'BackGroundColor',satellite_color)
    
    set(h4a,'MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);
    set(h4b,'MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);
    set(h4c,'MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);
    
elseif length(satellite_color) == 1
    
    satellite_color = get(handles.pushbutton3,'BackGroundColor');
    
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(~, ~, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Declare global variables
global h3a h3b h3c orbit_color

orbit_color     = uisetcolor;

if length(orbit_color) ~= 1
    
    set(handles.pushbutton4,'BackGroundColor',orbit_color)
    
    set(h3a,'Color',orbit_color);
    set(h3b,'Color',orbit_color);
    set(h3c,'Color',orbit_color);
    
elseif length(orbit_color) == 1
    
    orbit_color = get(handles.pushbutton4,'BackGroundColor');
    
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(~, ~, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Declare global variables
global h5a h5b velocity_color

velocity_color     = uisetcolor;

if length(velocity_color) ~= 1
    
    set(handles.pushbutton5,'BackGroundColor',velocity_color)
    
    set(h5a,'Color',velocity_color);
    set(h5b,'Color',velocity_color);
    
elseif length(velocity_color) == 1
    
    velocity_color = get(handles.pushbutton5,'BackGroundColor');
    
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(~, ~, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Declare global variables
global h6a h6b radial_color

radial_color     = uisetcolor;

if length(radial_color) ~= 1
    
    set(handles.pushbutton6,'BackGroundColor',radial_color)
    
    set(h6a,'Color',radial_color);
    set(h6b,'Color',radial_color);
    
elseif length(radial_color) == 1
    
    radial_color = get(handles.pushbutton6,'BackGroundColor');
    
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(~, ~, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Declare global variables
global h7a h7b irf_color

irf_color                                       = uisetcolor;

if length(irf_color) ~= 1
    
    set(handles.pushbutton7,'BackGroundColor',irf_color)
    
    CData_irf                                   = ones(2,21,3);
    CData_irf(:,:,1)                            = CData_irf(:,:,1).*irf_color(1);
    CData_irf(:,:,2)                            = CData_irf(:,:,2).*irf_color(2);
    CData_irf(:,:,3)                            = CData_irf(:,:,3).*irf_color(3);
    
    set(h7a(1),'Color',irf_color)
    set(h7a(2:4),'CData',CData_irf)
    set(h7b,'Color',irf_color);
    
elseif length(irf_color) == 1
    
    irf_color                                   = get(handles.pushbutton7,'BackGroundColor');
    
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(~, ~, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Declare global variables
global h8a h8b efrf_color

efrf_color     = uisetcolor;

if length(efrf_color) ~= 1
    
    set(handles.pushbutton8,'BackGroundColor',efrf_color)
    
    set(h8a,'Color',efrf_color);
    set(h8b,'Color',efrf_color);
    
elseif length(efrf_color) == 1
    
    efrf_color = get(handles.pushbutton8,'BackGroundColor');
    
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(~, ~, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UTC_offset = str2double(get(handles.edit13,'String'));

if UTC_offset >= 12
        
    errordlg('Please select a value between -12 and 12 hours.')
    
else
    
    set(handles.edit13,'String',num2str(UTC_offset + 0.25));
    show_current_time(0, 0, handles)
    
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(~, ~, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UTC_offset = str2double(get(handles.edit13,'String'));

if UTC_offset <= -12
        
    errordlg('Please select a value between -12 and 12 hours.')
    
else
    
    set(handles.edit13,'String',num2str(UTC_offset - 0.25));
    show_current_time(0, 0, handles)
    
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.checkbox15,'Value') == 1 || get(handles.checkbox16,'Value') == 1 || get(handles.checkbox17,'Value') == 1
    
    if get(handles.checkbox12,'Value') == 0
        
        animate_simulated_orbit(hObject, eventdata, handles)
        
    else
        
        animate_real_time_orbit(hObject, eventdata, handles)
                
    end
    
else
    
    warndlg('Please check at least one reference frame.')
    
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(~, ~, ~)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

satlab_DeleteFcn
satlab


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

checkbox11_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(~, ~, ~)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(~, ~, ~)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(~, ~, ~)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Menu_1_Callback(~, ~, ~)
% hObject    handle to Menu_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Menu_1_1_Callback(~, ~, handles)
% hObject    handle to Menu_1_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Declare global variables
global handles2

handles2 = handles;

catalog_selection


% --------------------------------------------------------------------
function Menu_1_2_Callback(~, ~, ~)
% hObject    handle to Menu_1_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

satellite_info


function show_current_time(~, ~, handles)

%Declare global variables
global current_epoch

c             = clock;
UTC_offset    = str2double(get(handles.edit13,'String'));

current_epoch = f_gd2jd(c(3),c(2),c(1),c(4)+UTC_offset,c(5),c(6),'JD');

set(handles.edit11,'String',num2str(current_epoch));

theta_z       = f_ut12lmst(0,current_epoch,'JD');

set(handles.edit12,'String',theta_z);


function orbit_real_time(~, ~, handles)

%Declare global variables
global POS_EARTH POS_COAST POS_t row col row2 col2 irf_axes efrf_axes
global GM period_Earth omega_Earth O a e w W i M theta_z T_num T t era
global h1a h2a h4a h5a h6a h8a h9a h3b h4b h5b h6b h7b h9b h3c h4c h9c
global sat_data_new sat_num current_epoch LATITUDE_t LONGITUDE_t day_night_map
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat phi_sat phi_sat_corr lamda_sat lamda_sat_corr
global X_EARTH_NEW Y_EARTH_NEW Z_EARTH_NEW X_COAST_NEW Y_COAST_NEW Z_COAST_NEW EFRF_AXES IRF_AXES X_t_NEW Y_t_NEW Z_t_NEW
global x_sat_c y_sat_c z_sat_c v_x_sat_c v_y_sat_c v_z_sat_c X_sat_c Y_sat_c Z_sat_c V_X_sat_c V_Y_sat_c V_Z_sat_c phi_sat_c lamda_sat_c
global f1_rt f2_rt f3_rt h1a_rt h2a_rt h4a_rt h5a_rt h6a_rt h8a_rt h9a_rt h3b_rt h4b_rt h5b_rt h6b_rt h7b_rt h9b_rt h3c_rt h4c_rt h9c_rt

%Get reference epoch in JD
t_0                                                     = f_gd2jd(sat_data_new(sat_num).epoch_day,1,sat_data_new(sat_num).epoch_year,0,0,0,'JD');

%Get current epoch in JD
t_c                                                     = current_epoch;

%Get current UTC epoch in day of year format
c                                                       = clock;

t_c_DOY                                                 = current_epoch - f_gd2jd(0,1,c(1),0,0,0,'JD');

%Calculate day night map
day_night_map                                           = f_day_night_map(t_c_DOY,LATITUDE_t,LONGITUDE_t);

%Propagate the orbit
theta_z                                                 = f_ut12lmst(0,t_c,'JD')*pi/180;
T_num                                                   = 1;
T                                                       = 2*pi*sqrt(a^3/GM)/period_Earth;
t                                                       = (t_c - (T/2 + T/360):T/360: t_c + T/2 + T/360)';
era                                                     = theta_z + omega_Earth*(t - t_c)*86400;

%Inertial Reference Frame
[x_sat_c,y_sat_c,z_sat_c,v_x_sat_c,v_y_sat_c,v_z_sat_c] = f_kepl2svec(GM,a,e,w,W,i,M,t_0,t_c);
[x_sat,y_sat,z_sat,v_x_sat,v_y_sat,v_z_sat]             = f_kepl2svec(GM,a,e,w,W,i,M,t_0,t);

%Initialize variables
X_sat                                                   = zeros(size(x_sat));
Y_sat                                                   = zeros(size(y_sat));
Z_sat                                                   = zeros(size(z_sat));

V_X_sat                                                 = zeros(size(v_x_sat));
V_Y_sat                                                 = zeros(size(v_y_sat));
V_Z_sat                                                 = zeros(size(v_z_sat));

%Earth Fixed Reference Frame
pos_vec_c                                               = [x_sat_c;y_sat_c;z_sat_c];
vel_vec_c                                               = [v_x_sat_c;v_y_sat_c;v_z_sat_c];

POS_VEC_c                                               = f_rotz(theta_z)*pos_vec_c;
VEL_VEC_c                                               = f_rotz(theta_z)*vel_vec_c - O*POS_VEC_c;

X_sat_c                                                 = POS_VEC_c(1,1);
Y_sat_c                                                 = POS_VEC_c(2,1);
Z_sat_c                                                 = POS_VEC_c(3,1);

V_X_sat_c                                               = VEL_VEC_c(1,1);
V_Y_sat_c                                               = VEL_VEC_c(2,1);
V_Z_sat_c                                               = VEL_VEC_c(3,1);

for ii = 1:size(x_sat,1)
    
    pos_vec                                             = [x_sat(ii,1);y_sat(ii,1);z_sat(ii,1)];
    vel_vec                                             = [v_x_sat(ii,1);v_y_sat(ii,1);v_z_sat(ii,1)];
    
    POS_VEC                                             = f_rotz(era(ii,1))*pos_vec;
    VEL_VEC                                             = f_rotz(era(ii,1))*vel_vec - O*POS_VEC;
    
    X_sat(ii,1)                                         = POS_VEC(1,1);
    Y_sat(ii,1)                                         = POS_VEC(2,1);
    Z_sat(ii,1)                                         = POS_VEC(3,1);
    
    V_X_sat(ii,1)                                       = VEL_VEC(1,1);
    V_Y_sat(ii,1)                                       = VEL_VEC(2,1);
    V_Z_sat(ii,1)                                       = VEL_VEC(3,1);
    
end

POS_EARTH_NEW                                           = f_rotz(-theta_z)*POS_EARTH;
POS_COAST_NEW                                           = f_rotz(-theta_z)*POS_COAST;
POS_t_NEW                                               = f_rotz(-theta_z)*POS_t;

X_EARTH_NEW                                             = reshape(POS_EARTH_NEW(1,:),row,col);
Y_EARTH_NEW                                             = reshape(POS_EARTH_NEW(2,:),row,col);
Z_EARTH_NEW                                             = reshape(POS_EARTH_NEW(3,:),row,col);

X_COAST_NEW                                             = POS_COAST_NEW(1,:)';
Y_COAST_NEW                                             = POS_COAST_NEW(2,:)';
Z_COAST_NEW                                             = POS_COAST_NEW(3,:)';

X_t_NEW                                                 = reshape(POS_t_NEW(1,:),row2,col2);
Y_t_NEW                                                 = reshape(POS_t_NEW(2,:),row2,col2);
Z_t_NEW                                                 = reshape(POS_t_NEW(3,:),row2,col2);

IRF_AXES                                                = f_rotz(-theta_z)*irf_axes;
EFRF_AXES                                               = f_rotz(theta_z)*efrf_axes;

%Earth Fixed Reference Frame (Map Projection)
[phi_sat_c,lamda_sat_c,~]                               = f_cart2geo(X_sat_c,Y_sat_c,Z_sat_c,'GRS 1980');
[phi_sat,lamda_sat,~]                                   = f_cart2geo(X_sat,Y_sat,Z_sat,'GRS 1980');

%Fix east-west lines bug
[phi_sat_corr,lamda_sat_corr]                           = f_fix_lines(phi_sat,lamda_sat,200);

%Plot to Inertial Reference Frame
set(h1a,'XData',X_EARTH_NEW,'YData',Y_EARTH_NEW,'ZData',Z_EARTH_NEW);
set(h2a,'XData',X_COAST_NEW,'YData',Y_COAST_NEW,'ZData',Z_COAST_NEW);
set(h4a,'XData',x_sat_c,'YData',y_sat_c,'ZData',z_sat_c);
set(h5a,'XData',x_sat_c,'YData',y_sat_c,'ZData',z_sat_c,'UData',v_x_sat_c*1000,'VData',v_y_sat_c*1000,'WData',v_z_sat_c*1000);
set(h6a,'XData',[0 x_sat_c],'YData',[0 y_sat_c],'ZData',[0 z_sat_c]);
set(h8a,'XData',[0;0;0],'YData',[0;0;0],'ZData',[0;0;0],'UData',EFRF_AXES(:,1),'VData',EFRF_AXES(:,2),'WData',EFRF_AXES(:,3));
set(h9a,'XData',X_t_NEW,'YData',Y_t_NEW,'ZData',Z_t_NEW,'CData',day_night_map);
set(h9a,'Visible','on')

%Plot to Earth Fixed Reference Frame
set(h3b,'XData',X_sat,'YData',Y_sat,'ZData',Z_sat);
set(h4b,'XData',X_sat_c,'YData',Y_sat_c,'ZData',Z_sat_c);
set(h5b,'XData',X_sat_c,'YData',Y_sat_c,'ZData',Z_sat_c,'UData',V_X_sat_c*1000,'VData',V_Y_sat_c*1000,'WData',V_Z_sat_c*1000);
set(h6b,'XData',[0 X_sat_c],'YData',[0 Y_sat_c],'ZData',[0 Z_sat_c]);
set(h7b,'XData',[0;0;0],'YData',[0;0;0],'ZData',[0;0;0],'UData',IRF_AXES(:,1),'VData',IRF_AXES(:,2),'WData',IRF_AXES(:,3));
set(h9b,'CData',day_night_map);
set(h9b,'Visible','on')

%Plot to Earth Fixed Reference Frame (Map Projection)
set(h3c,'XData',lamda_sat_corr,'YData',phi_sat_corr);
set(h4c,'XData',lamda_sat_c,'YData',phi_sat_c);
set(h9c,'XData',LONGITUDE_t,'YData',LATITUDE_t,'CData',day_night_map);
set(h9c,'Visible','on')

%Plot to figures

%Plot to Inertial Reference Frame
if get(handles.checkbox15,'Value') == 1 & ishandle(f1_rt) == 1
    
    set(h1a_rt,'XData',X_EARTH_NEW,'YData',Y_EARTH_NEW,'ZData',Z_EARTH_NEW);
    set(h2a_rt,'XData',X_COAST_NEW,'YData',Y_COAST_NEW,'ZData',Z_COAST_NEW);
    set(h4a_rt,'XData',x_sat_c,'YData',y_sat_c,'ZData',z_sat_c);
    set(h5a_rt,'XData',x_sat_c,'YData',y_sat_c,'ZData',z_sat_c,'UData',v_x_sat_c*1000,'VData',v_y_sat_c*1000,'WData',v_z_sat_c*1000);
    set(h6a_rt,'XData',[0 x_sat_c],'YData',[0 y_sat_c],'ZData',[0 z_sat_c]);
    set(h8a_rt,'XData',[0;0;0],'YData',[0;0;0],'ZData',[0;0;0],'UData',EFRF_AXES(:,1),'VData',EFRF_AXES(:,2),'WData',EFRF_AXES(:,3));
    set(h9a_rt,'XData',X_t_NEW,'YData',Y_t_NEW,'ZData',Z_t_NEW,'CData',day_night_map);
    
end

%Plot to Earth Fixed Reference Frame
if get(handles.checkbox16,'Value') == 1 & ishandle(f2_rt) == 1
    
    set(h3b_rt,'XData',X_sat,'YData',Y_sat,'ZData',Z_sat);
    set(h4b_rt,'XData',X_sat_c,'YData',Y_sat_c,'ZData',Z_sat_c);
    set(h5b_rt,'XData',X_sat_c,'YData',Y_sat_c,'ZData',Z_sat_c,'UData',V_X_sat_c*1000,'VData',V_Y_sat_c*1000,'WData',V_Z_sat_c*1000);
    set(h6b_rt,'XData',[0 X_sat_c],'YData',[0 Y_sat_c],'ZData',[0 Z_sat_c]);
    set(h7b_rt,'XData',[0;0;0],'YData',[0;0;0],'ZData',[0;0;0],'UData',IRF_AXES(:,1),'VData',IRF_AXES(:,2),'WData',IRF_AXES(:,3));
    set(h9b_rt,'CData',day_night_map);
    
end

%Plot to Earth Fixed Reference Frame (Map Projection)
if get(handles.checkbox17,'Value') == 1 & get(handles.checkbox5,'Value') == 1 & ishandle(f3_rt) == 1
    
    set(h3c_rt,'XData',lamda_sat_corr,'YData',phi_sat_corr);
    set(h4c_rt,'XData',lamda_sat_c,'YData',phi_sat_c);
    set(h9c_rt,'XData',LONGITUDE_t,'YData',LATITUDE_t,'CData',day_night_map);
    
end


function animate_simulated_orbit(~, ~, handles)

%Declare global variables
global POS_EARTH POS_COAST LON LAT TOPOGRAPHY row col long_corr lat_corr irf_axes efrf_axes
global theta_z t era
global earth_color coastline_color satellite_color orbit_color velocity_color radial_color irf_color efrf_color
global x_sat y_sat z_sat v_x_sat v_y_sat v_z_sat X_sat Y_sat Z_sat V_X_sat V_Y_sat V_Z_sat phi_sat phi_sat_corr lamda_sat lamda_sat_corr

%Calculate initial states
X_EARTH                                     = reshape(POS_EARTH(1,:),row,col);
Y_EARTH                                     = reshape(POS_EARTH(2,:),row,col);
Z_EARTH                                     = reshape(POS_EARTH(3,:),row,col);

X_COAST                                     = POS_COAST(1,:)';
Y_COAST                                     = POS_COAST(2,:)';
Z_COAST                                     = POS_COAST(3,:)';

POS_EARTH_NEW                               = f_rotz(-theta_z)*POS_EARTH;
POS_COAST_NEW                               = f_rotz(-theta_z)*POS_COAST;

X_EARTH_NEW                                 = reshape(POS_EARTH_NEW(1,:),row,col);
Y_EARTH_NEW                                 = reshape(POS_EARTH_NEW(2,:),row,col);
Z_EARTH_NEW                                 = reshape(POS_EARTH_NEW(3,:),row,col);

X_COAST_NEW                                 = POS_COAST_NEW(1,:)';
Y_COAST_NEW                                 = POS_COAST_NEW(2,:)';
Z_COAST_NEW                                 = POS_COAST_NEW(3,:)';

IRF_AXES                                    = f_rotz(-theta_z)*irf_axes;
EFRF_AXES                                   = f_rotz(theta_z)*efrf_axes;

%Find subplot number
subplot_num                                 = 0;
subplot_cur                                 = 0;

if get(handles.checkbox15,'Value') == 1
    
    subplot_num                             = subplot_num + 1;
    
end

if get(handles.checkbox16,'Value') == 1
    
    subplot_num                             = subplot_num + 1;
    
end

if get(handles.checkbox17,'Value') == 1
    
    subplot_num                             = subplot_num + 1;
    
end

%Get screen and figure dimensions
screen_dimensions                           = get(0,'screensize');
screen_width                                = screen_dimensions(3);
screen_height                               = screen_dimensions(4);
figure_width                                = screen_width/subplot_num;
figure_height                               = screen_height;

%Plot to Inertial Reference Frame
if get(handles.checkbox15,'Value') == 1
    
    subplot_cur                             = subplot_cur + 1;
    
    f1                                      = figure;
    
    set(f1,'Position',[subplot_cur*(screen_width/subplot_num) - figure_width, 100, figure_width, figure_height - 200]);
    
    %Plot Earth
    h1a                                     = surf(X_EARTH_NEW,Y_EARTH_NEW,Z_EARTH_NEW,TOPOGRAPHY);
    
    if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
        
        demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
        shading interp
        
    elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
        
        set(h1a,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color);
        
    elseif get(handles.checkbox2,'Value') == 0
        
        set(h1a,'Visible','off')
        
    end
    
    hold on
    
    %Plot coastline
    h2a                                     = plot3(X_COAST_NEW,Y_COAST_NEW,Z_COAST_NEW,'Color',coastline_color);
    
    if get(handles.checkbox4,'Value') == 0
        
        set(h2a,'Visible','off')
        
    end
    
    %Plot satellite orbit
    h3a                                     = plot3(x_sat(1:362),y_sat(1:362),z_sat(1:362),'Color',orbit_color,'Linewidth',2);
    
    if get(handles.checkbox6,'Value') == 0
        
        set(h3a,'Visible','off')
        
    end
    
    %Plot satellite
    h4a                                     = scatter3(x_sat(1,1),y_sat(1,1),z_sat(1,1),100,'filled','MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);
    
    if get(handles.checkbox5,'Value') == 0
        
        set(h4a,'Visible','off')
        
    end
    
    %Plot satellite velocity
    h5a                                     = quiver3(x_sat(1,1),y_sat(1,1),z_sat(1,1),v_x_sat(1,1)*1000,v_y_sat(1,1)*1000,v_z_sat(1,1)*1000,'Color',velocity_color,'Linewidth',2);
    
    if get(handles.checkbox7,'Value') == 0
        
        set(h5a,'Visible','off')
        
    end
    
    %Plot satellite radial distance
    h6a                                     = plot3([0 x_sat(1,1)],[0 y_sat(1,1)],[0 z_sat(1,1)],'Color',radial_color);
    
    if get(handles.checkbox8,'Value') == 0
        
        set(h6a,'Visible','off')
        
    end
    
    %Plot Inertial Reference Frame axes
    h7a                                     = quiver3([0;0;0],[0;0;0],[0;0;0],irf_axes(:,1),irf_axes(:,2),irf_axes(:,3),'Color',irf_color,'Linewidth',2,'ShowArrowHead','off');
    
    if get(handles.checkbox9,'Value') == 0
        
        set(h7a,'Visible','off')
        
    end
    
    %Plot Earth Fixed Reference Frame axes
    h8a                                     = quiver3([0;0;0],[0;0;0],[0;0;0],EFRF_AXES(:,1),EFRF_AXES(:,2),EFRF_AXES(:,3),'Color',efrf_color,'Linewidth',2,'ShowArrowHead','off');
    
    if get(handles.checkbox10,'Value') == 0
        
        set(h8a,'Visible','off')
        
    end
    
    axis equal
    axis off
    grid off
    box off
    rotate3d on
    hold off
    
end

%Plot to Earth Fixed Reference Frame
if get(handles.checkbox16,'Value') == 1
    
    subplot_cur                             = subplot_cur + 1;
    
    f2                                      = figure;
    
    set(f2,'Position',[subplot_cur*(screen_width/subplot_num) - figure_width, 100, figure_width, figure_height - 200]);
    
    %Plot Earth
    h1b                                     = surf(X_EARTH,Y_EARTH,Z_EARTH,TOPOGRAPHY);
    
    if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
        
        demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
        shading interp
        
    elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
        
        set(h1b,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color);
        
    elseif get(handles.checkbox2,'Value') == 0
        
        set(h1b,'Visible','off')
        
    end
    
    hold on
    
    %Plot coastline
    h2b                                     = plot3(X_COAST,Y_COAST,Z_COAST,'Color',coastline_color);
    
    if get(handles.checkbox4,'Value') == 0
        
        set(h2b,'Visible','off')
        
    end
    
    %Plot satellite orbit
    h3b                                     = plot3(X_sat,Y_sat,Z_sat,'Color',orbit_color,'Linewidth',2);
    
    if get(handles.checkbox6,'Value') == 0
        
        set(h3b,'Visible','off')
        
    end
    
    %Plot satellite
    h4b                                     = scatter3(X_sat(1,1),Y_sat(1,1),Z_sat(1,1),100,'filled','MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);
    
    if get(handles.checkbox5,'Value') == 0
        
        set(h4b,'Visible','off')
        
    end
    
    %Plot satellite velocity
    h5b                                     = quiver3(X_sat(1,1),Y_sat(1,1),Z_sat(1,1),V_X_sat(1,1)*1000,V_Y_sat(1,1)*1000,V_Z_sat(1,1)*1000,'Color',velocity_color,'Linewidth',2);
    
    if get(handles.checkbox7,'Value') == 0
        
        set(h5b,'Visible','off')
        
    end
    
    %Plot satellite radial distance
    h6b                                     = plot3([0 X_sat(1,1)],[0 Y_sat(1,1)],[0 Z_sat(1,1)],'Color',radial_color);
    
    if get(handles.checkbox8,'Value') == 0
        
        set(h6b,'Visible','off')
        
    end
    
    %Plot Inertial Reference Frame axes
    h7b                                     = quiver3([0;0;0],[0;0;0],[0;0;0],IRF_AXES(:,1),IRF_AXES(:,2),IRF_AXES(:,3),'Color',irf_color,'Linewidth',2,'ShowArrowHead','off');
    
    if get(handles.checkbox9,'Value') == 0
        
        set(h7b,'Visible','off')
        
    end
    
    %Plot Earth Fixed Reference Frame axes
    h8b                                     = quiver3([0;0;0],[0;0;0],[0;0;0],efrf_axes(:,1),efrf_axes(:,2),efrf_axes(:,3),'Color',efrf_color,'Linewidth',2,'ShowArrowHead','off');
    
    if get(handles.checkbox10,'Value') == 0
        
        set(h8b,'Visible','off')
        
    end
    
    view(-24,48)
    axis equal
    axis off
    grid off
    box off
    rotate3d on
    hold off
    
end

%Plot to Earth Fixed Reference Frame (Ground Tracks)
if get(handles.checkbox17,'Value') == 1
    
    subplot_cur                             = subplot_cur + 1;
    
    f3                                      = figure;
    
    set(f3,'Position',[subplot_cur*(screen_width/subplot_num) - figure_width, 100, figure_width, figure_height - 200]);
    
    %Plot Earth
    h1c                                     = pcolor(LON,LAT,TOPOGRAPHY);
    
    if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
        
        demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
        shading interp
        
    elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
        
        set(h1c,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color);
        
    elseif get(handles.checkbox2,'Value') == 0
        
        set(h1c,'Visible','off')
        
    end
    
    hold on
    
    %Plot coastline
    h2c                                     = plot(long_corr,lat_corr,'Color',coastline_color);
    
    if get(handles.checkbox4,'Value') == 0
        
        set(h2c,'Visible','off')
        
    end
    
    %Plot satellite orbit
    h3c                                     = plot(lamda_sat_corr,phi_sat_corr,'Color',orbit_color,'Linewidth',2);
    
    if get(handles.checkbox6,'Value') == 0
        
        set(h3c,'Visible','off')
        
    end
    
    %Plot satellite
    h4c                                     = scatter(lamda_sat(1,1),phi_sat(1,1),100,'filled','MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);
    
    if get(handles.checkbox5,'Value') == 0
        
        set(h4c,'Visible','off')
        
    end
    
    axis equal
    axis off
    grid off
    box off
    hold off
    
end

%Update the plot data
for ii = 2:size(t,1)
    
    %Calculate states for each epoch
    POS_EARTH_NEW                           = f_rotz(-era(ii,1))*POS_EARTH;
    POS_COAST_NEW                           = f_rotz(-era(ii,1))*POS_COAST;
    
    X_EARTH_NEW                             = reshape(POS_EARTH_NEW(1,:),row,col);
    Y_EARTH_NEW                             = reshape(POS_EARTH_NEW(2,:),row,col);
    Z_EARTH_NEW                             = reshape(POS_EARTH_NEW(3,:),row,col);
    
    X_COAST_NEW                             = POS_COAST_NEW(1,:)';
    Y_COAST_NEW                             = POS_COAST_NEW(2,:)';
    Z_COAST_NEW                             = POS_COAST_NEW(3,:)';
    
    IRF_AXES                                = f_rotz(-era(ii,1))*irf_axes;
    EFRF_AXES                               = f_rotz(era(ii,1))*efrf_axes;
    
    %Plot to Inertial Reference Frame
    if get(handles.checkbox15,'Value') == 1 && ishandle(f1) == 1
        
        set(h1a,'XData',X_EARTH_NEW,'YData',Y_EARTH_NEW,'ZData',Z_EARTH_NEW);
        set(h2a,'XData',X_COAST_NEW,'YData',Y_COAST_NEW,'ZData',Z_COAST_NEW);
        
        if get(handles.checkbox5,'Value') == 1
            
            set(h4a,'XData',x_sat(ii,:),'YData',y_sat(ii,:),'ZData',z_sat(ii,:));
            
        end
        
        set(h5a,'XData',x_sat(ii,:),'YData',y_sat(ii,:),'ZData',z_sat(ii,:),'UData',v_x_sat(ii,:)*1000,'VData',v_y_sat(ii,:)*1000,'WData',v_z_sat(ii,:)*1000);
        set(h6a,'XData',[0,x_sat(ii,:)],'YData',[0,y_sat(ii,:)],'ZData',[0,z_sat(ii,:)]);
        set(h8a,'XData',[0;0;0],'YData',[0;0;0],'ZData',[0;0;0],'UData',EFRF_AXES(:,1),'VData',EFRF_AXES(:,2),'WData',EFRF_AXES(:,3));
        
    end
    
    %Plot to Earth Fixed Reference Frame
    if get(handles.checkbox16,'Value') == 1 && ishandle(f2) == 1
        
        if get(handles.checkbox5,'Value') == 1
            
            set(h4b,'XData',X_sat(ii,1),'YData',Y_sat(ii,1),'ZData',Z_sat(ii,1));
            
        end
        
        set(h5b,'XData',X_sat(ii,1),'YData',Y_sat(ii,1),'ZData',Z_sat(ii,1),'UData',V_X_sat(ii,1)*1000,'VData',V_Y_sat(ii,1)*1000,'WData',V_Z_sat(ii,1)*1000);
        set(h6b,'XData',[0 X_sat(ii,1)],'YData',[0 Y_sat(ii,1)],'ZData',[0 Z_sat(ii,1)]);
        set(h7b,'XData',[0;0;0],'YData',[0;0;0],'ZData',[0;0;0],'UData',IRF_AXES(:,1),'VData',IRF_AXES(:,2),'WData',IRF_AXES(:,3));
        
    end
    
    %Plot to Earth Fixed Reference Frame (Ground Tracks)
    if get(handles.checkbox17,'Value') == 1 && get(handles.checkbox5,'Value') == 1 && ishandle(f3) == 1
        
        set(h4c,'XData',lamda_sat(ii,1),'YData',phi_sat(ii,1),'ZData',1);
        
    end
    
    drawnow
    
end


function animate_real_time_orbit(~, ~, handles)

%Declare global variables
global irf_axes efrf_axes
global LATITUDE_t LONGITUDE_t day_night_map
global x_sat y_sat z_sat X_sat Y_sat Z_sat phi_sat_corr lamda_sat_corr
global long_corr lat_corr TOPOGRAPHY X_EARTH Y_EARTH Z_EARTH X_COAST Y_COAST Z_COAST X_t Y_t Z_t
global LON LAT X_EARTH_NEW Y_EARTH_NEW Z_EARTH_NEW X_COAST_NEW Y_COAST_NEW Z_COAST_NEW EFRF_AXES IRF_AXES X_t_NEW Y_t_NEW Z_t_NEW
global x_sat_c y_sat_c z_sat_c v_x_sat_c v_y_sat_c v_z_sat_c X_sat_c Y_sat_c Z_sat_c V_X_sat_c V_Y_sat_c V_Z_sat_c phi_sat_c lamda_sat_c
global earth_color coastline_color satellite_color orbit_color velocity_color radial_color irf_color efrf_color
global f1_rt f2_rt f3_rt h1a_rt h2a_rt h3a_rt h4a_rt h5a_rt h6a_rt h7a_rt h8a_rt h9a_rt h1b_rt h2b_rt h3b_rt h4b_rt h5b_rt h6b_rt h7b_rt h8b_rt h9b_rt h1c_rt h2c_rt h3c_rt h4c_rt h9c_rt

%Find subplot number
subplot_num                                 = 0;
subplot_cur                                 = 0;

if get(handles.checkbox15,'Value') == 1
    
    subplot_num                             = subplot_num + 1;
    
end

if get(handles.checkbox16,'Value') == 1
    
    subplot_num                             = subplot_num + 1;
    
end

if get(handles.checkbox17,'Value') == 1
    
    subplot_num                             = subplot_num + 1;
    
end

%Get screen and figure dimensions
screen_dimensions                           = get(0,'screensize');
screen_width                                = screen_dimensions(3);
screen_height                               = screen_dimensions(4);
figure_width                                = screen_width/subplot_num;
figure_height                               = screen_height;

%Plot to Inertial Reference Frame
if get(handles.checkbox15,'Value') == 1
    
    subplot_cur                             = subplot_cur + 1;
    
    f1_rt                                   = figure;
    
    set(f1_rt,'Position',[subplot_cur*(screen_width/subplot_num) - figure_width, 100, figure_width, figure_height - 200]);
    
    %Plot Earth
    h1a_rt                                  = surf(X_EARTH_NEW,Y_EARTH_NEW,Z_EARTH_NEW,TOPOGRAPHY);
    
    if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
        
        demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
        shading interp
        
    elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
        
        set(h1a_rt,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color);
        
    elseif get(handles.checkbox2,'Value') == 0
        
        set(h1a_rt,'Visible','off')
        
    end
    
    hold on
    
    %Plot day night map
    h9a_rt                                  = surf(X_t_NEW,Y_t_NEW,Z_t_NEW,day_night_map);
    
    shading interp
    alpha(h9a_rt,0.5)
    
    %Plot coastline
    h2a_rt                                  = plot3(X_COAST_NEW,Y_COAST_NEW,Z_COAST_NEW,'Color',coastline_color);
    
    if get(handles.checkbox4,'Value') == 0
        
        set(h2a_rt,'Visible','off')
        
    end
    
    %Plot satellite orbit
    h3a_rt                                  = plot3(x_sat(1:362),y_sat(1:362),z_sat(1:362),'Color',orbit_color,'Linewidth',2);
    
    if get(handles.checkbox6,'Value') == 0
        
        set(h3a_rt,'Visible','off')
        
    end
    
    %Plot satellite
    h4a_rt                                  = scatter3(x_sat_c,y_sat_c,z_sat_c,100,'filled','MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);
    
    if get(handles.checkbox5,'Value') == 0
        
        set(h4a_rt,'Visible','off')
        
    end
    
    %Plot satellite velocity
    h5a_rt                                  = quiver3(x_sat_c,y_sat_c,z_sat_c,v_x_sat_c*1000,v_y_sat_c*1000,v_z_sat_c*1000,'Color',velocity_color,'Linewidth',2);
    
    if get(handles.checkbox7,'Value') == 0
        
        set(h5a_rt,'Visible','off')
        
    end
    
    %Plot satellite radial distance
    h6a_rt                                  = plot3([0 x_sat_c],[0 y_sat_c],[0 z_sat_c],'Color',radial_color);
    
    if get(handles.checkbox8,'Value') == 0
        
        set(h6a_rt,'Visible','off')
        
    end
    
    %Plot Inertial Reference Frame axes
    h7a_rt                                  = quiver3([0;0;0],[0;0;0],[0;0;0],irf_axes(:,1),irf_axes(:,2),irf_axes(:,3),'Color',irf_color,'Linewidth',2,'ShowArrowHead','off');
    
    if get(handles.checkbox9,'Value') == 0
        
        set(h7a_rt,'Visible','off')
        
    end
    
    %Plot Earth Fixed Reference Frame axes
    h8a_rt                                  = quiver3([0;0;0],[0;0;0],[0;0;0],EFRF_AXES(:,1),EFRF_AXES(:,2),EFRF_AXES(:,3),'Color',efrf_color,'Linewidth',2,'ShowArrowHead','off');
    
    if get(handles.checkbox10,'Value') == 0
        
        set(h8a_rt,'Visible','off')
        
    end
    
    axis equal
    axis off
    grid off
    box off
    rotate3d on
    hold off
    
end

%Plot to Earth Fixed Reference Frame
if get(handles.checkbox16,'Value') == 1
    
    subplot_cur                             = subplot_cur + 1;
    
    f2_rt                                   = figure;
    
    set(f2_rt,'Position',[subplot_cur*(screen_width/subplot_num) - figure_width, 100, figure_width, figure_height - 200]);
    
    %Plot Earth
    h1b_rt                                  = surf(X_EARTH,Y_EARTH,Z_EARTH,TOPOGRAPHY);
    
    if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
        
        demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
        shading interp
        
    elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
        
        set(h1b_rt,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color);
        
    elseif get(handles.checkbox2,'Value') == 0
        
        set(h1b_rt,'Visible','off')
        
    end
    
    hold on
    
    %Plot day night map
    h9b_rt                                  = surf(X_t,Y_t,Z_t,day_night_map);
    
    shading interp
    alpha(h9b_rt,0.5)
    
    %Plot coastline
    h2b_rt                                  = plot3(X_COAST,Y_COAST,Z_COAST,'Color',coastline_color);
    
    if get(handles.checkbox4,'Value') == 0
        
        set(h2b_rt,'Visible','off')
        
    end
    
    %Plot satellite orbit
    h3b_rt                                  = plot3(X_sat,Y_sat,Z_sat,'Color',orbit_color,'Linewidth',2);
    
    if get(handles.checkbox6,'Value') == 0
        
        set(h3b_rt,'Visible','off')
        
    end
    
    %Plot satellite
    h4b_rt                                  = scatter3(X_sat_c,Y_sat_c,Z_sat_c,100,'filled','MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);
    
    if get(handles.checkbox5,'Value') == 0
        
        set(h4b_rt,'Visible','off')
        
    end
    
    %Plot satellite velocity
    h5b_rt                                  = quiver3(X_sat_c,Y_sat_c,Z_sat_c,V_X_sat_c*1000,V_Y_sat_c*1000,V_Z_sat_c*1000,'Color',velocity_color,'Linewidth',2);
    
    if get(handles.checkbox7,'Value') == 0
        
        set(h5b_rt,'Visible','off')
        
    end
    
    %Plot satellite radial distance
    h6b_rt                                  = plot3([0 X_sat_c],[0 Y_sat_c],[0 Z_sat_c],'Color',radial_color);
    
    if get(handles.checkbox8,'Value') == 0
        
        set(h6b_rt,'Visible','off')
        
    end
    
    %Plot Inertial Reference Frame axes
    h7b_rt                                  = quiver3([0;0;0],[0;0;0],[0;0;0],IRF_AXES(:,1),IRF_AXES(:,2),IRF_AXES(:,3),'Color',irf_color,'Linewidth',2,'ShowArrowHead','off');
    
    if get(handles.checkbox9,'Value') == 0
        
        set(h7b_rt,'Visible','off')
        
    end
    
    %Plot Earth Fixed Reference Frame axes
    h8b_rt                                  = quiver3([0;0;0],[0;0;0],[0;0;0],efrf_axes(:,1),efrf_axes(:,2),efrf_axes(:,3),'Color',efrf_color,'Linewidth',2,'ShowArrowHead','off');
    
    if get(handles.checkbox10,'Value') == 0
        
        set(h8b_rt,'Visible','off')
        
    end
    
    axis equal
    axis off
    grid off
    box off
    rotate3d on
    hold off
    
end

%Plot to Earth Fixed Reference Frame (Ground Tracks)
if get(handles.checkbox17,'Value') == 1
    
    subplot_cur                             = subplot_cur + 1;
    
    f3_rt                                   = figure;
    
    set(f3_rt,'Position',[subplot_cur*(screen_width/subplot_num) - figure_width, 100, figure_width, figure_height - 200]);
    
    %Plot Earth
    h1c_rt                                  = pcolor(LON,LAT,TOPOGRAPHY);
    
    if get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 1
        
        demcmap([min(min(TOPOGRAPHY)) max(max(TOPOGRAPHY))]);
        shading interp
        
    elseif get(handles.checkbox2,'Value') == 1 && get(handles.checkbox3,'Value') == 0
        
        set(h1c_rt,'CData',ones(size(TOPOGRAPHY)),'EdgeColor','none','FaceColor',earth_color);
        
    elseif get(handles.checkbox2,'Value') == 0
        
        set(h1c_rt,'Visible','off')
        
    end
    
    hold on
    
    %Plot day night map
    h9c_rt                                  = pcolor(LONGITUDE_t,LATITUDE_t,day_night_map);
    
    shading interp
    alpha(h9c_rt,0.5)
    
    %Plot coastline
    h2c_rt                                  = plot(long_corr,lat_corr,'Color',coastline_color);
    
    if get(handles.checkbox4,'Value') == 0
        
        set(h2c_rt,'Visible','off')
        
    end
    
    %Plot satellite orbit
    h3c_rt                                  = plot(lamda_sat_corr,phi_sat_corr,'Color',orbit_color,'Linewidth',2);
    
    if get(handles.checkbox6,'Value') == 0
        
        set(h3c_rt,'Visible','off')
        
    end
    
    %Plot satellite
    h4c_rt                                  = scatter(lamda_sat_c,phi_sat_c,100,'filled','MarkerEdgeColor',satellite_color,'MarkerFaceColor',satellite_color);
    
    set(h4c_rt,'ZData',10)
    
    if get(handles.checkbox5,'Value') == 0
        
        set(h4c_rt,'Visible','off')
        
    end
    
    axis equal
    axis off
    grid off
    box off
    hold off
    
end
