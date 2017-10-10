load('data.mat','B1','B3','B4');
B = {B1,B2,B3,B4};
t = toepoch([2001 9 10 7 50 0]);
DATABASE = c_ctl(0,'isdat_db');
% B1 = B5vps;
% c_eval('[~,~,B1]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');',jc);
% % B1dsi=c_coord_trans('GSE','DSC',B1,'cl_id',jc);
thetacell = {0,0,0,0};
phicell = {0,0,0,0};
for jc = [1 2 3 4]
ts = linspace(t,t+480,length(B{jc}))';
%%   Get phase data, phase_data is in degrees.
% [time phase] = caa_is_get(DATABASE,t-5,485,1,'ephemeris','phase_2');
% phase = phase*pi/180; % to radians
% phase = [time phase]; % make array
% phase = c_phase(ts,phase); % interpolate
% time = phase(:,1); % define new time vector

%% Get phase data from c_load
[ok,phase]=c_load('Atwo?',jc);
phase(:,2) = phase(:,2)*pi/180;
phase = c_phase(ts,phase); % interpolates to ts resolution
time = phase(:,1);

%%   Convert B to DSI frame
Bdsi = c_coord_trans('GSE','DSI',B{jc},'cl_id',jc);

%%   Normalization of B
% take only x and y components:
bxy = Bdsi(:,1:3);

bn = irf_norm(bxy(:,2:3)); % must specify that only columns 2 and 3 should be normalized
bn3 = irf_norm(Bdsi(:,2:4));
bx = bn(:,1);
by = bn(:,2);

bx3 = bn3(:,1);
by3 = bn3(:,2);
bz3 = bn3(:,3);

%%   Convert from phase angle to direction vector for P12.
p12 = [cos(phase(:,2)-pi/4) -sin(phase(:,2)-pi/4)]; 
p12 = [p12 zeros(length(p12),1)];
%% Make time series
theta = acos(bx.*p12(:,1) + by.*p12(:,2))*180/pi;
theta = [time theta];

%% Make absolute angle
phi = acos(bx3.*p12(:,1) + by3.*p12(:,2) + bz3.*p12(:,3))*180/pi;
phi = [time phi];
% make cell array for use by spectest.m
thetacell{jc} = theta;
phicell{jc} = phi;
%%   Plot the time series
% phase(:,2) = phase(:,2)*180/pi;
% figure
% h1(1) = irf_subplot(2,1,-1);
% irf_plot(phase);
% ylabel('Phase [deg]')
% titlestring = sprintf('Phase and B-field angle for Cluster %d',jc);
% title(titlestring);
% h1(2) = irf_subplot(2,1,-2);
% irf_plot(theta);
% irf_zoom(h1,'x',[t t+60]);
% add_timeaxis(h1(2),'date');
% ylabel('Angle B-p12 [deg]');
end