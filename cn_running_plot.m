if 0 % Loading data
c_eval('[caaB?,~,gseB?fgm]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');',3:4);
%% Magnetic field
c_eval('gsmB?=irf_abs(gsmB?);',3:4);
c_eval('gseB?=c_coord_trans(''GSM'',''GSE'',gsmB?,''cl_id'',?);',3:4);
c_eval('gseB?=irf_abs(gseB?);',3:4);
%% Electric field
c_eval('[caaE?,~,diE?]=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'');',3:4);
c_eval('gseE?=c_coord_trans(''dsi'',''gse'',diE?,''CL_ID'',?);',3:4);
%% Spacecraft potential
c_eval('[caaP?,~,P?]=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'');',3:4);
%% Electron densities
c_eval('scpNe?=c_efw_scp2ne(P?);',3:4);
c_eval('[caaPNe?,~,peaNe?]=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'');',3:4);
c_eval('peaNe?=irf_resamp(peaNe?,diE?);',3:4);
%% Ion densities
c_eval('[caaNi?,~,Ni?]=c_caa_var_get(''density__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',3);
%% Ion and electron (ExB) velocities GSE
c_eval('[caaExB?,~,diExB?]=c_caa_var_get(''v_drift_ISR2__C?_CP_EFW_L2_V3D_INERT'');',3:4);
c_eval('gseExB?=c_coord_trans(''dsi'',''gse'',diExB?,''CL_ID'',?);',3:4);
c_eval('[caaVi?,~,hiaVi?]=c_caa_var_get(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',3);
c_eval('[caaVi?,~,codifVi?]=c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'');',3:4 );
c_eval('hiaVi?=irf_resamp(hiaVi?,diE3);',3);
c_eval('codifVi?=irf_resamp(codifVi?,diE3);',3:4);
%% Spacecraft positions
c_eval('[caaPos?,~,gsePos?]=c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_FULL'');',3:4);         
%% Electron temperature
% K 10^6
c_eval('[caaparTe?,~,parTe?]=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'');',3:4);
c_eval('[caaperTe?,~,perTe?]=c_caa_var_get(''Data_Temperature_ComponentPerpendicularToMagField__C?_CP_PEA_'');',3:4);
% Adding to one temperature only, using C3
Te3=[parTe3(:,1) (parTe3(:,2)+2*perTe3(:,2))/3];


% eV
eVTe3=[Te3(:,1) Te3(:,2)*8.61734*10];
eVTe3=irf_resamp(eVTe3,diE3);

%Te4=[parTe4(:,1) (parTe4(:,2)+2*perTe4(:,2))/3];
%eVTe4=[Te4(:,1) Te4(:,2)*8.61734*10];
%eVTe4=irf_resamp(eVTe4,diE3);
%Tetot4=[Tepar4(:,1) (Tepar4(:,2)+2*Teper4(:,2))/3];
%Teper=cn_average(Teper3,Teper4);
%% Ion temperatures
% K 10^6
[caaTi4,~,Ti4]=c_caa_var_get('temperature__C3_CP_CIS_HIA_ONBOARD_MOMENTS');
% eV
eVTi4=[Ti4(:,1) Ti4(:,2)*8.61734e-2*1000];
eVTi4=irf_resamp(eVTi4,diE3);
end
% 0831: CODIF finns p? b?da, HIA p? C3
% 0902: samma
%% Time intervals
if 1 
%% 20070831 10:18:40.60 - 2007 08 31 10 10:18:42.00
t1    =[2007 08 31 10 18 40 60]; 
t2    =[2007 08 31 10 18 42 00];
n_hat=-[0.0856    0.9293   -0.3594];

t1    =[2007 08 31 10 18 41 00]; 
t2    =[2007 08 31 10 18 42 00];
n_hat=-[0.0856    0.9293   -0.3594];

t1    =[2007 08 31 10 18 41 20]; 
t2    =[2007 08 31 10 18 41 50];
n_hat=-[0.0856    0.9293   -0.3594];
%% 20070831 10:18:42.50 - 2007 08 31 10 10:18:43.00
t1    =[2007 08 31 10 18 42 50];
t2    =[2007 08 31 10 18 43 00];
n_hat=-[0.0856    0.9293   -0.3594];
%% 20070831 10:18:46.00 - 2007 08 31 10 10:18:46.50
t1    =[2007 08 31 10 18 45 80];
t2    =[2007 08 31 10 18 46 60];
n_hat=-[0.0856    0.9293   -0.3594];
%% 20070831 10:18:46.00 - 2007 08 31 10 10:18:46.50
t1    =[2007 08 31 10 19 02 95];
t2    =[2007 08 31 10 19 03 40];
n_hat=-[0.0856    0.9293   -0.3594];
%% 20070831 10:19:02.60 - 2007 08 31 10 10:19:03.50
t1    =[2007 08 31 10 19 02 60];
t2    =[2007 08 31 10 19 03 40];
n_hat=-[0.0856    0.9293   -0.3594];
%n_hat=[-0.1470 -0.9490 0.2789];
%% 20070831 10:19:02.60 - 2007 08 31 10 10:19:03.00 X
t1    =[2007 08 31 10 19 02.60];
t2    =[2007 08 31 10 19 03.00];
n_hat=-[0.0856    0.9293   -0.3594];

%n_hat=[-0.1470 -0.9490 0.2789];
%% 20070831 10:19:03.00 - 2007 08 31 10 10:19:03.40
t1    =[2007 08 31 10 19 03 00];
t2    =[2007 08 31 10 19 03 40];
n_hat=[0.0856    0.9293   -0.3594];

n_hat=[-0.1470 -0.9490 0.2789]; % stor dubbelpuls
%% 20070831 10:19:03.40 - 2007 08 31 10 10:19:04.20
t1    =[2007 08 31 10 19 03.50];
t2    =[2007 08 31 10 19 04.50];
n_hat=[0.0856    0.9293   -0.3594];
%% 20070831 10:19:03.40 - 2007 08 31 10 10:19:04.20
t1    =[2007 08 31 10 19 03 40];
t2    =[2007 08 31 10 19 04 20];
n_hat=[0.0856    0.9293   -0.3594];
%% 20070831 10:19:04.50 - 2007 08 31 10 10:19:05.10 X
t1    =[2007 08 31 10 19 04.50]; t2    =[2007 08 31 10 19 05.10];
ta    =[2007 08 31 10 19 04.70]; tb    =[2007 08 31 10 19 04.90];
n_hat=-[0.0856    0.9293   -0.3594]; % global value
%n_hat=[0.1280   -0.9860    0.1066]; % local value
%n_hat=-[0.2315    0.8372   -0.4955]; % from cn_tool
%n_hat=[ 0.065 -0.973  0.221];
%% 20070831 10:19:05.50 - 2007 08 31 10 10:19:05.90 X
t1    =[2007 08 31 10 19 05.50];
t2    =[2007 08 31 10 19 05.90];
n_hat=-[0.0856    0.9293   -0.3594];
%% 20070831 10:19:05.95 - 2007 08 31 10 10:19:06.55
t1    =[2007 08 31 10 19 05.95];
t2    =[2007 08 31 10 19 06.55];
n_hat=[0.0856    0.9293   -0.3594];
%% 20070831 10:19:05.95 - 2007 08 31 10 10:19:06.55 X
t1    =[2007 08 31 10 19 06.15];
t2    =[2007 08 31 10 19 06.40];
n_hat=-[0.0856    0.9293   -0.3594];
%% 20070831 10:19:06.50 - 2007 08 31 10 10:19:08.00 X not in report though 
t1    =[2007 08 31 10 19 06.90];
t2    =[2007 08 31 10 19 07.50];
n_hat=-[0.0856    0.9293   -0.3594];
%%
t1    =[2007 08 31 10 19 06.90];
t2    =[2007 08 31 10 19 07.25];
n_hat=-[0.0856    0.9293   -0.3594];
%%
t1    =[2007 08 31 10 19 07.20];
t2    =[2007 08 31 10 19 07.50];
n_hat=-[0.0856    0.9293   -0.3594];
%% 20070831 10:19:10.12 - 2007 08 31 10 10:19:11.10 X not in report though
t1    =[2007 08 31 10 19 10.12];
t2    =[2007 08 31 10 19 11.10];
n_hat=-[0.0856    0.9293   -0.3594];
%%
t1    =[2007 08 31 10 19 10.12];
t2    =[2007 08 31 10 19 10.45];
n_hat=-[0.0856    0.9293   -0.3594];
%%
t1    =[2007 08 31 10 19 10.60];
t2    =[2007 08 31 10 19 11.10];
n_hat=-[0.0856    0.9293   -0.3594];

%% 20070831 10:39:28.00 - 2007 08 31 10 10:39:30.00 too low frequency...
t1    =[2007 08 31 10 39 28 00];
t2    =[2007 08 31 10 39 30 00];
n_hat=[-0.0001    0.9501    0.3120];
n_hat=[-0.0778   -0.9543   -0.2886];
n_hat=[-0.3578   -0.9086   -0.2156]; % only on very local area
%%
if 0
%% 2007 08 31 10 19 03 30
tstart=[2007 08 31 10 19 01 00];
tstop =[2007 08 31 10 19 05 00];
t1    =[2007 08 31 10 19 03 00];
t2    =[2007 08 31 10 19 03 30];

n_hat=[0.0856    0.9293   -0.3594]; % from bigger area
%% 2007 08 31 10 19 07 30
tstart=[2007 08 31 10 19 06 50];
tstop =[2007 08 31 10 19 08 00];
t1    =[2007 08 31 10 19 07 50];
t2    =[2007 08 31 10 19 08 00];
n_hat=[0.0856    0.9293   -0.3594];
%% 2007 08 31 10 19 03 30 DL?
tstart=[2007 08 31 10 19 08 00];
tstop =[2007 08 31 10 19 05 00];
t1    =[2007 08 31 10 19 10 50];
t2    =[2007 08 31 10 19 11 50];

n_hat=[0.0297 0.9459 -0.3231];
%% 2007 08 31 10 18 40 05
tstart=[2007 08 31 10 18 44 10];
tstop =[2007 08 31 10 18 48 20];
t1    =[2007 08 31 10 18 45 10];
t2    =[2007 08 31 10 18 45 20];

n_hat=[0.0856    0.9293   -0.3594]; % from bigger area
%% 2007 08 31 10 18 42 45
tstart=[2007 08 31 10 18 41 00];
tstop =[2007 08 31 10 18 45 00];
t1    =[2007 08 31 10 18 42 45];
t2    =[2007 08 31 10 18 43 55];

n_hat=[0.2989 0.9128 0.2783]; % current layer normal in GSE
n_hat=[0.2552 0.9275 0.2732];
n_hat=-[-0.3919 0.8750 -0.2841];
n_hat=[0 1 0]; % ExB-direction
n_hat=[0.0910 0.9260 -0.3665];
end

%% %%%%%%%%%%%%%%%%%%%% 2 september %%%%%%%%%%%%%%%%%%%%%%%%%%%


t1    =[2007 09 02 17 40 00 00];
t2    =[2007 09 02 18 10 00 00];
%%
t1    =[2007 09 02 15 20 00 00];
t2    =[2007 09 02 15 35 30 00];
%%
t1    =[2007 09 02 15 45 00 00];
t2    =[2007 09 02 15 50 30 00];
%% 20070902 15:47:00.00 - 2007 09 02 15:49:00.00
t1    =[2007 09 02 15 30 20 00];
t2    =[2007 09 02 15 30 30 00];
n_hat=[0.1255    0.3195   -0.9392]; %L2/L3=6.6 Bmin=0
% inte lika p? C3 och C4
%% 20070902 15:47:00.00 - 2007 09 02 15:49:00.00
t1    =[2007 09 02 15 33 10 00];
t2    =[2007 09 02 15 33 30 00];



n_hat=[-0.2704    0.9544   -0.1267]; %L2/L3=7 Bmin=2
% inte lika p? C3 och C4
%% 20070902 15:47:00.00 - 2007 09 02 15:49:00.00
t1    =[2007 09 02 15 47 00 00];
t2    =[2007 09 02 15 48 30 00];
n_hat=[0.0856    0.9293   -0.3594]; % fel
%% 20070902 15:47:14.60 - 2007 09 02 15:47:16.20 ok
t1    =[2007 09 02 15 47 14 00];
t2    =[2007 09 02 15 47 16 00];
n_hat=[0.3064    0.5762   -0.7577];
%% 20070902 15:47:00.00 - 2007 09 02 15:49:00.00
t1    =[2007 09 02 15 47 17 30];
t2    =[2007 09 02 15 47 18 50];
n_hat=[-0.3496   -0.5600    0.7511]; %L2/L3=87 Bmin=1
%% 20070902 15:47:31.60 - 2007 09 02 15:47:33.20 ok
t1    =[2007 09 02 15 47 31 60];
t2    =[2007 09 02 15 47 33 50];
n_hat=[-0.5399   0.4423    0.7162]; % L2/L3=8
%% 20070902 15:47:31.60 - 2007 09 02 15:47:33.20 X
t1    =[2007 09 02 15 47 29 30];
t2    =[2007 09 02 15 47 30 80];
n_hat=[-0.5704    0.4338    0.6975]; %L3/L2=37 Bmin=4
%n_hat=[-0.590  0.710  0.385];
%%
t1    =[2007 09 02 14 32 29 00];
t2    =[2007 09 02 14 32 30 00];
%%
t1    =[2007 09 02 14 32 29 00];
t2    =[2007 09 02 14 32 32 00];
%%
t1    =[2007 09 02 14 32 31 00];
t2    =[2007 09 02 14 32 32 00];
%%
t1    =[2007 09 02 14 32 35 00];
t2    =[2007 09 02 14 32 36 00];
%%
t1    =[2007 09 02 14 34 07 00];
t2    =[2007 09 02 14 34 08 00];
%% E finns
t1    =[2007 09 02 14 34 54 60];
t2    =[2007 09 02 14 34 55 10];
%% E finns
t1    =[2007 09 02 14 34 55 10];
t2    =[2007 09 02 14 34 55 80];
%%
t1    =[2007 09 02 14 32 29 00];
t2    =[2007 09 02 14 32 30 00];
%% 20070902 15:47:31.60 - 2007 09 02 15:47:33.20 ok X
t1    =[2007 09 02 15 47 32 00];
t2    =[2007 09 02 15 47 33 40];
n_hat=[-0.5399   0.4423    0.7162]; % L2/L3=8 Bmin=0
%% 20070902 15:48:04.00 - 2007 09 02 15:48:12.00
t1    =[2007 09 02 15 48 05 00];
t2    =[2007 09 02 15 48 90 00];
n_hat=[-0.2237   -0.7907    0.5699]; %L2/L3=74 Bmin=0
%% 20070902 15:48:04.00 - 2007 09 02 15:48:12.00
t1    =[2007 09 02 15 48 05 50];
t2    =[2007 09 02 15 48 07 50];
n_hat=[0.2443    0.8458   -0.4743];
%% 20070902 15:48:04.00 - 2007 09 02 15:48:12.00
t1    =[2007 09 02 15 47 29 00];
t2    =[2007 09 02 15 47 31 00];
n_hat=[-0.4976    0.4924    0.7141];
%% 20070902 15:48:04.00 - 2007 09 02 15:48:12.00
t1    =[2007 09 02 15 47 31 10];
t2    =[2007 09 02 15 47 34 00];
n_hat=[-0.5564    0.4338    0.7087];
%% 20070902 15:48:04.00 - 2007 09 02 15:48:12.00
t1    =[2007 09 02 17 58 16 10];
t2    =[2007 09 02 17 58 18 00];
n_hat=[-0.3929    0.7131   -0.5806]; % L2/L3=98 Bmin =2


%% %%%%%%%%%%%%%% 26 september %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1=[2007 09 26 10 40 00 00];
t2=[2007 09 26 11 00 00 00];

%% 20070926 15:48:04.00 - 2007 09 26 15:48:12.00
t1    =[2007 09 26 10 50 25 20];
t2    =[2007 09 26 10 50 25 90];
n_hat=-[0.2443    0.8458   -0.4743];
%n_hat=[0.6770    0.0158   -0.7359]; % L2/L3=7.5 Bmin=3
%% 20070926 15:48:04.00 - 2007 09 26 15:48:12.00
t1    =[2007 09 26 10 49 10 00];
t2    =[2007 09 26 10 49 25 00];
n_hat=[0.3231   -0.0598   -0.9445];
%% 20070926 15:48:04.00 - 2007 09 26 15:48:12.00
t1    =[2007 09 26 10 17 24 00];
t2    =[2007 09 26 10 17 28 00];
n_hat=[0.1079   -0.6849   -0.7206];
%% 20070926 10:51:17.00 - 2007 09 26 10:51:18.20
t1    =[2007 09 26 10 51 17 00];
t2    =[2007 09 26 10 51 18 20];
n_hat=[-0.6629   -0.6999    0.2659];
%% 20070926 10:51:17.00 - 2007 09 26 10:51:18.20
t1    =[2007 09 26 10 50 08 00];
t2    =[2007 09 26 10 50 12 20];
n_hat=[0.4971   -0.2348   -0.8353];
%% 20070926 10:51:17.00 - 2007 09 26 10:51:18.20
t1    =[2007 09 26 10 50 25 20];
t2    =[2007 09 26 10 50 26 00];
n_hat=[0.4971   -0.2348   -0.8353];
%% 20070926 10:51:17.00 - 2007 09 26 10:51:18.20
t1    =[2007 09 26 10 50 27 80];
t2    =[2007 09 26 10 50 28 20];
n_hat=[-0.8007    0.1327    0.5842];
%% 20070926 10:51:17.00 - 2007 09 26 10:51:18.20
t1    =[2007 09 26 10 16 34 00];
t2    =[2007 09 26 10 16 34 80];
n_hat=[0.0719   -0.4893   -0.8692]; % max_var EAC-direction=vd
%% 20070926 10:51:17.00 - 2007 09 26 10:51:18.20
t1    =[2007 09 26 09 49 07 30];
t2    =[2007 09 26 09 49 08 00];
n_hat=[0.0719   -0.4893   -0.8692];

%%  
t1    =[2007 09 26 10 45 00 00];
t2    =[2007 09 26 10 55 00 00];
end
%%
t1str=[num2str(t1(1)),'0',num2str(t1(2)),'',num2str(t1(3)),...
            '',num2str(t1(4)),'',num2str(t1(5)),'',num2str(t1(6))];
t2str=[num2str(t2(1)),'0',num2str(t2(2)),'',num2str(t2(3)),...
            '',num2str(t2(4)),'',num2str(t2(5)),'',num2str(t2(6))];

% Taking average density during time interval
c_eval('peaNe?z=cn_toepoch(t1,t2,peaNe?);',3:4);
c_eval('peaNe?av=sum(peaNe?z(:,2),1)/size(peaNe?z,1);',3:4);
peaNeav=(peaNe3av+peaNe4av)/2;

Niav=cn_toepoch(t1,hiaNi3);

%c_eval('Ne?z=cn_toepoch(t1,t2,Ne?);',3:4)
%c_eval('Ne?z=[Ne?z(:,1) Ne?z(:,2)/10];',3:4)

d=cn_event_conf(t1,t2,gseB3,gseB4,codifVi3,codifVi4,gseExB3,gseExB4,gsePos3,gsePos4,eVTe3,eVTi4,peaNeav,Niav,scpNe3,peaNe3,n_hat,t1str,t2str,'poster');
%%
% Transform into usable quantities when using nhat like basis vector
Bhat=d{1};
M=d{2};
Vi=d{3};
vihat=d{4};
B=d{5};
vdhat=d{6}; % possible propagation direction of wave in field aligned coord
vd=(M\vdhat')'; % possible propagation direction of wave in GSE coord
bh=Bhat;
flh=d{7};
Vith=d{8};
Tiav=d{9};
Teav=d{10};
normdist=d{11};
kdist=d{12};
r_i=d{13};
zdist=d{14};
r_e=d{15};
Ln1=d{16};
% Check1 existance condition again, there is another funtion in
% cn_event_conf
Ln2=cn_ln2(t1,t2,t1str,t2str,gsmB3,gsmB4,Teav,Tiav,peaNeav,Niav,normdist,kdist,M,r_i);
%%
if 1 % Reducing time series
    tone=t1;
    ttwo=t2;
    B3zoom=cn_toepoch(tone,ttwo,gseB3);
    B4zoom=cn_toepoch(tone,ttwo,gseB4);
    E3zoom=cn_toepoch(tone,ttwo,diE3);
    E4zoom=cn_toepoch(tone,ttwo,diE4);
    ExB3zoom=cn_toepoch(tone,ttwo,gseExB3);
    ExB4zoom=cn_toepoch(tone,ttwo,gseExB4);
    c_eval('ExB?abs=irf_abs(ExB?zoom);',3:4);
    % Getting them to be the same length 
    if size(E3zoom,1)>size(E4zoom,1); E4zoom=irf_resamp(E4zoom,E3zoom);
    elseif size(E3zoom,1)<size(E4zoom,1); E3zoom=irf_resamp(E3zoom,E4zoom);
    end
end
if 1 % Separating DC and AC E-field
dt=E3zoom(2,1)-E3zoom(1,1);
fs=1/dt;
flow=flh*0.5;
fhigh=180;
% Taking a slightly larger interval to make the filtering
t11=t1;
t11(5)=t11(5)-1; % subtracting one minute
t22=t2;
t22(5)=t22(5)+1; % adding one minute

c_eval('E?zoom2=cn_toepoch(t11,t22,diE?);',3:4)

if size(E3zoom2,1)>size(E4zoom2,1); E4zoom2=irf_resamp(E4zoom2,E3zoom2);
elseif size(E3zoom2,1)<size(E4zoom2,1); E3zoom2=irf_resamp(E3zoom2,E4zoom2);
end

c_eval('E?AC=irf_filt(E?zoom(:,1:3),flow,0,fs,3);',3:4);
c_eval('E?DC=irf_filt(E?zoom(:,1:3),0,flow,fs,3);',3:4);
c_eval('E?DC2=[E?zoom(:,1) E?zoom(:,[2 3])-E?AC(:,[2 3])];',3:4);

c_eval('E?AC=irf_filt(E?zoom2(:,1:3),flow,0,fs,3);',3:4);
c_eval('E?DC=irf_filt(E?zoom2(:,1:3),0,flow,fs,3);',3:4);
c_eval('E?DC2=[E?zoom2(:,1) E?zoom2(:,[2 3])-E?AC(:,[2 3])];',3:4);
c_eval('E?AC2=[E?zoom2(:,1) E?zoom2(:,[2 3])-E?DC(:,[2 3])];',3:4)

c_eval('E?AC=cn_toepoch(t1,t2,E?AC);',3:4)
c_eval('E?AC2=cn_toepoch(t1,t2,E?AC2);',3:4)
c_eval('E?DC=cn_toepoch(t1,t2,E?DC);',3:4)
c_eval('E?DC2=cn_toepoch(t1,t2,E?DC2);',3:4)

end 
if 0 % Plot E (DCAC DC AC)
figure('name','E: tot DC AC')
h=irf_plot(10);
isub=1;
if 1 % Allt
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3zoom(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4zoom(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);  
        irf_zoom(hca,'y');       

        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3zoom(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4zoom(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{Y,}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);   
        irf_zoom(hca,'y'); 
end
if 1 % DC direkt
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3DC(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4DC(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X,DC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);       

        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3DC(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4DC(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{Y,DC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);  
end  
if 1 % Allt minus AC=DC   
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3DC2(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4DC2(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X,DC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);   
        irf_zoom(hca,'y');      

        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3DC2(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4DC2(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{Y,DC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);   
        irf_zoom(hca,'y'); 
end  
if 1 % AC direkt
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3AC(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4AC(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X,AC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);  
        irf_zoom(hca,'y');  
        
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3AC(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4AC(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{Y,AC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);   
        irf_zoom(hca,'y');
end    
if 1 % Allt minus DC=AC
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3AC2(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4AC2(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X,AC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);  
        irf_zoom(hca,'y');  
        
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3AC(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4AC(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{Y,AC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);   
        irf_zoom(hca,'y');
end  

    title(h(1),['DC/AC E-field (ISR2)      cut at ',num2str(flow),'Hz'])
    irf_zoom(h,'x',[cn_toepoch(t1) cn_toepoch(t2)])
    irf_plot_axis_align
    irf_timeaxis(hca,'usefig');
    set(h,'YLimMode','auto');
    %
    eval(['print -dpdf ',t1str,'_',t2str,'_E_ACDC.pdf']);
end

% The one I will us I overwrite as AC? and DC?
c_eval('E?AC=E?AC2;',3:4)
c_eval('E?DC=E?DC;',3:4)

%
f=1; % x/y/xy
g=3; % DCAC/DC/AC
window=10;
if 1 %%%% Getting correlation of E-field using mycorr
corrx=zeros(size(E3zoom,1)*2-1,3);
corry=zeros(size(E3zoom,1)*2-1,3);
corrxy=zeros(size(E3zoom,1)*2-1,3);

isub=1;
[t0x(isub) dtx(isub) tplusx(isub) tminusx(isub) corrx(:,isub)]=cn_mycorr(E3zoom(:,2),E4zoom(:,2),window,fs);
[t0y(isub) dty(isub) tplusy(isub) tminusy(isub) corry(:,isub)]=cn_mycorr(E3zoom(:,3),E4zoom(:,3),window,fs);
[t0x(isub) dtxy(isub) tplusxy(isub) tminusxy(isub) corrxy(:,isub)]=cn_mycorr(E3zoom(:,2:3),E4zoom(:,2:3),window,fs); 
isub=isub+1;
% Getting correlation for DC-field
[t0x(isub) dtx(isub) tplusx(isub) tminusx(isub) corrx(:,isub)]=cn_mycorr(E3DC(:,2),E4DC(:,2),window,fs); 
[t0y(isub) dty(isub) tplusy(isub) tminusy(isub) corry(:,isub)]=cn_mycorr(E3DC(:,3),E4DC(:,3),window,fs);
[t0xy(isub) dtxy(isub) tplusxy(isub) tminusxy(isub) corrxy(:,isub)]=cn_mycorr(E3DC(:,2:3),E4DC(:,2:3),window,fs);
isub=isub+1;
% Getting correlation for AC-field
[t0x(isub) dtx(isub) tplusx(isub) tminusx(isub) corrx(:,isub)]=cn_mycorr(E3AC(:,2),E4AC(:,2),window,fs);
[t0y(isub) dty(isub) tplusy(isub) tminusy(isub) corry(:,isub)]=cn_mycorr(E3AC(:,3),E4AC(:,3),window,fs);
[t0xy(isub) dtxy(isub) tplusxy(isub) tminusxy(isub) corrxy(:,isub)]=cn_mycorr(E3AC(:,2:3),E4AC(:,2:3),window,fs);
end
if 0 % plot correlation
 figure
 hca=subplot(3,1,1);
 plot(hca,corrx(:,1),'b');hold on;   
 plot(hca,corry(:,1),'g'); hold on;   
 plot(hca,corrxy(:,1),'c'); hold on;  
 legend('x','y','x+y')
 ylabel(hca,'DC+AC')
 hca=subplot(3,1,2);
 title('Correlation for C3 and C4 E-field  (mycorr)')      
 plot(hca,corrx(:,2),'b');hold on;   
 plot(hca,corry(:,2),'g'); hold on;   
 plot(hca,corrxy(:,2),'c'); hold on;  
 legend(hca,'x','y','x+y')
 ylabel(hca,'DC')   
 hca=subplot(3,1,isub);   
 plot(hca,corrx(:,3),'b');hold on;
 plot(hca,corry(:,3),'g'); hold on;
 plot(hca,corrxy(:,3),'c'); hold on;
 isub=isub+1;
 legend('x','y','x+y')
 ylabel(hca,'AC')
 eval(['print -dpdf ',t1str,'_',t2str,'_mycorr.pdf']);
end
if 1 % Getting velocity from t0 (mycorr)
t0matrix=[t0x;t0y;t0xy]; %
dtmatrix=[dtx;dty;dtxy];
for j=1:3 % x, y, xy
    for k=1:3 % DC+AC, DC, AC
        t0=t0matrix(j,k);
        dt=dtmatrix(j,k);
        [v(j,k) vminus(j,k) vplus(j,k)]=cn_v(t0,dt,vd,gsePos3,gsePos4,(t1+t2)/2);
    end
end
end
if 1 % Correlate field in propagation direction
        
end
if 1 %%%% Using crosscorrelation matlab function %%%%%%%%%%%%%%%%%%%%%%
corrx2=zeros(size(E3zoom,1)*2-1,3);
corry2=zeros(size(E3zoom,1)*2-1,3);
corrxy2=zeros(size(E3zoom,1)*2-1,3);

isub=1;
[t0x2(isub) corrx2(:,isub)]=cn_xcorr(E3zoom,E4zoom,window,'x');
[t0y2(isub) corry2(:,isub)]=cn_xcorr(E3zoom,E4zoom,window,'y');
[t0xy2(isub) corrxy2(:,isub)]=cn_xcorr(E3zoom,E4zoom,window,'xy');
isub=isub+1;
% Getting correlation for DC-field
[t0x2(isub) corrx2(:,isub)]=cn_xcorr(E3DC,E4DC,window,'x');
[t0y2(isub) corry2(:,isub)]=cn_xcorr(E3DC,E4DC,window,'y');
[t0xy2(isub) corrxy2(:,isub)]=cn_xcorr(E3DC,E4DC,window,'xy');
isub=isub+1;
% Getting correlation for AC-field
[t0x2(isub) corrx2(:,isub)]=cn_xcorr(E3AC,E4AC,window,'x');
[t0y2(isub) corry2(:,isub)]=cn_xcorr(E3AC,E4AC,window,'y');
[t0xy2(isub) corrxy2(:,isub)]=cn_xcorr(E3AC,E4AC,window,'xy');
end
if 0 % plot correlation
figure
 hca=subplot(3,1,1);
 plot(hca,corrx2(:,1),'b');hold on;
 plot(hca,corry2(:,1),'g'); hold on;
 plot(hca,corrxy2(:,1),'c'); hold on;
 legend('x','y','x+y')
 ylabel(hca,'DC+AC')
 title('Correlation for C3 and C4 E-field (xcorr)') 
 hca=subplot(3,1,2);
 plot(hca,corrx2(:,2),'b');hold on;
 plot(hca,corry2(:,2),'g'); hold on;
 plot(hca,corrxy2(:,2),'c'); hold on;
 legend(hca,'x','y','x+y')
 ylabel(hca,'DC')
 hca=subplot(3,1,3);
 plot(hca,corrx2(:,3),'b');hold on;
 plot(hca,corry2(:,3),'g'); hold on;
 plot(hca,corrxy2(:,3),'c'); hold on;  
 isub=isub+1;
 legend('x','y','x+y')
 ylabel(hca,'AC')
 eval(['print -dpdf ',t1str,'_',t2str,'_xcorr.pdf']);
end
if 1 % Getting velocity from t0 (xcorr)
clear t0matrix2 v2;
t0matrix2=[t0x2; t0y2; t0xy2];
dtmatrix2=[t0x2; t0y2; t0xy2];

for j=1:3  % x, y, xy
    for k=1:3 % DC+AC, DC, AC
        t02=t0matrix2(j,k);
        dt2=dtmatrix2(j,k);
        [v2(j,k),~,~]=cn_v(t02,dt2,vd,gsePos3,gsePos4,(t1+t2)/2);
    end
end
end
t0choice=t0matrix(f,g);
vchoice=v(f,g);
if 0
%% Getting velocity automatically
t0=t0x(3);
dt=dtx(3);
[v vminus vplus]=cn_v(t0,dt,vd,gsePos3,gsePos4,(t1+t2)/2)
end
if 0 % Plot shifted E field, for mycorr 
    cn_plot_corr(E3zoom,E4zoom,t0matrix(f,g),'mycorr');
    eval(['print -dpdf ',t1str,'_',t2str,'_E_mycorr.pdf']);
end
if 0 % and xcorr
    cn_plot_corr(E3zoom,E4zoom,t0matrix2(f,g),'xcorr');
    eval(['print -dpdf ',t1str,'_',t2str,'_E_xcorr.pdf']);
end


vchoice=vchoice;
t0choice=t0choice;
% Perpendicular plot
[bB3AC bB4AC]=cn_edir4(t1,t2,t1,t2,gseE3,gseE4,gseExB3,gseExB4,gsmB3,gsmB4,gsePos3,gsePos4,vchoice,M,flh,n_hat,Teav,t0choice,'AC',12);
%%
%t0choice=-t0choice;
%vchoice=vchoice;
if 1
    %c_eval('gseExB?=irf_resamp(gseExB?,gseE?);',3:4)
    c_eval('v_exb?=irf_filt(cn_toepoch(t1,t2,gseExB?),0,flow,fs,5);',3:4)
    c_eval('v_exb_n?=(M(1,:)*v_exb?(:,2:4)'')'';',3:4)
    c_eval('v_vec?=vchoice*repmat(vd,size(v_exb_n?,1),1);',3:4)
    c_eval('v_vec?(:,1)=v_vec?(:,1)-v_exb_n?;',3:4)
    c_eval('v_vec_t?=[v_exb?(:,1) v_vec?];',3:4)
    c_eval('gseB?AC=irf_filt(cn_toepoch(t1,t2,gseB?),flow,0,fs,5);',3:4)
    %gseB3AC=staffB3AC;
    %gseB4AC=staffB4AC;
    c_eval('bB?=cn_m_trans(gseB?AC,M,1);',3:4)
    c_eval('gseE?AC=irf_filt(cn_toepoch(t1,t2,gseE?),flow,0,fs,5);',3:4)
    %c_eval('phi?=cn_potential(cn_toepoch(t1,t2,gseE?AC),''GSE'',vd*vchoice,''GSE'');',3:4)
    c_eval('phi?=cn_potential(cn_toepoch(t1,t2,gseE?AC),''GSE'',v_vec_t?,''GSE'');',3:4)
    c_eval('bE?=cn_m_trans(cn_toepoch(t1,t2,gseE?),M,1);',3:4)
    c_eval('bE?AC=cn_m_trans(cn_toepoch(t1,t2,gseE?AC),M,1);',3:4)
    c_eval('gsmB?AC=irf_filt(cn_toepoch(t1,t2,gsmB?),flow,0,fs,5);',3:4)
    c_eval('gsmB?DC=irf_filt(cn_toepoch(t1,t2,gsmB?),0,flow,fs,5);',3:4)
    Ln3=cn_ln2(t1,t2,t1str,t2str,gsmB3DC,gsmB4DC,Teav,Tiav,peaNeav,Niav,normdist,kdist,M,r_i);
    % correlate e-field in propagation direction
    [t0k dtk tplusk tminusk corrk]=cn_mycorr(bE3AC(:,3),bE4AC(:,3),window,fs);
    t0choice=t0k;
    vchoice=kdist/t0k;
    vplus=kdist/tplusk;
    vminus=kdist/tminusk;
    %vpar=zdist/t0choice;
    %c_eval('phi?par=cn_potential(gseE?AC)',3:4)
end

%[t0k dtk tplusk tminusk corrk]=cn_mycorr(bE3AC(:,3),bE4AC(:,3),window,fs);

% Normalized potential
%[phiV3 phiV4]=cn_phi_plot(E3zoom,E4zoom,E3AC,E4AC,(t1+t2)/2,vd,vchoice,Teav,000,t1str,t2str,'AC',t0choice);
%[phi3 phi4]=cn_phi_plot2(E3zoom,E4zoom,E3AC,E4AC,(t1+t2)/2,vd,vchoice,Teav,t0choice);
%
%phi=phi3;
%gseB3AC=cn_toepoch(t1,t2,gseB3);
%bB3=cn_m_trans(gseB3AC,M,1);
% Compare dB and potential
Cstr='4';
comp_dB_phi_article

Cstr='3';
comp_dB_phi_article

propdir='rightleft';
propdir='topbot';

cn_ebplot2(bE3AC,bE4AC,bB3,bB4,...
    gsePos3,gsePos4,vchoice,M,propdir,r_e,n_hat,12);
% Parallel fields



%%
if 0
%% Calculating phi for whole E
c_eval('E?_toint=irf_abs(E?zoom);',3:4);
c_eval('phi?t=irf_integrate(E?_toint(:,1:3));',3:4);
c_eval('phi?=[phi?t(:,1) -phi?t(:,2:3).*vvec];',3:4);
%% Calculating phi for E, DC
c_eval('E?_toint=irf_abs(E?DC2);',3:4);
c_eval('phi?t=irf_integrate(E?_toint(:,1:3));',3:4);
c_eval('phi?=[phi?t(:,1) -phi?t(:,2:3).*vvec];',3:4);
end
if 0
%% Calculating phi for E, AC
c_eval('E?_toint=irf_abs(E?AC);',3:4);
c_eval('phi?t=irf_integrate(E?_toint(:,1:3));',3:4);
c_eval('phi?=[phi?t(:,1) -phi?t(:,2:3).*v?vec];',3:4);
end
% Plot potential with electric field and
% click on figure to obtain the potential drop
if 0
figure;

h=irf_plot(4);
isub=1;
if 1 % Electric field, AC
    hca=h(isub);isub=isub+1;
    irf_plot(hca,E3AC(:,[1 2]),'g'); hold(hca,'on');
    irf_plot(hca,E4AC(:,[1 2]),'b'); hold(hca,'on');
    ylabel(hca,'E_{X,DC}[mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);       

    hca=h(isub);isub=isub+1;
    irf_plot(hca,E3AC(:,[1 3]),'g'); hold(hca,'on');
    irf_plot(hca,E4AC(:,[1 3]),'b'); hold(hca,'on');
    ylabel(hca,'E_{Y,DC}[mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);  
end
if 1 % Potential
    hca=h(isub);isub=isub+1;
    irf_plot(hca,phi3(:,[1 2]),'g'); hold(hca,'on');
    irf_plot(hca,phi4(:,[1 2]),'b'); hold(hca,'on');
    ylabel(hca,'\phi_{X}[mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);       

    hca=h(isub);isub=isub+1;
    irf_plot(hca,phi3(:,[1 3]),'g'); hold(hca,'on');
    irf_plot(hca,phi4(:,[1 3]),'b'); hold(hca,'on');
    ylabel(hca,'\phi_{Y}[mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);  
end

title(h(1),'E-AC-field and electrostatic potential \phi (ISR2)')
irf_zoom(h,'x',[cn_toepoch(t1) cn_toepoch(t2)])
irf_plot_axis_align
add_timeaxis(hca,'usefig');

% Extracting time difference for potential drop 
 
[t pot]=ginput(4);
potx=pot(2)-pot(1)
poty=pot(4)-pot(3)

end
%% Saturation amplitudes

% Drake: ignored nonlinear coupling into damped modes with fintite k_par
% and is therefore likely to be an overestimate of A_sat

Vi_perp=sqrt(Vi(1)^2+Vi(2)^2);
phi_sat=0.0792*Vi_perp/Vith;

% Gary: assumed Vith>V which is not the case here
dE=35;


% Ion trapping V/Vith~>3 which is more the case here...
 

if 0
%%
vbl_gse=M\[vbl 0 0]'; % ion velocity due to boundary layer moving
vd_gse=M\[0 vd 0]'; % ion drift velocity
vparl_gse=M\[0 0 vpar]'; % field aligned ion velocity
%%
vph_hat=vd_gse/cn_mag(vd_gse); % propagation direction vector for LHD waves
                               % use in c_4_v_gui()
%%
time=toepoch([2007 08 31 10 18 43]);
vph_hat_dsi=c_coord_trans('dsi','gse',[time vph_hat'],'CL_ID',3);                               

%%
c_eval('gseE?=c_coord_trans(''dsi'',''gse'',[toepoch([2007 08 31 10 18 43]) [1 0 -1]],''CL_ID'',?);',3:4);
end