tint = irf.tint('2016-01-01T00:00:00.000Z/2017-01-01T00:00:00.000Z');
tint = irf.tint('2019-01-01T00:00:00.000Z/2020-01-01T00:00:00.000Z');
tint = irf.tint('2015-01-01T00:00:00.000Z/2016-01-01T00:00:00.000Z');


% K->eV conversion factor
units = irf_units; 
K2eV = units.kB/units.eV;

%pbm_orig = irf_get_data(tint,'n,n,P,B,Ms','omni_min');
%omni_data_vxyz = irf_get_data_omni(tint,'vx,vy,vz','omni2');

% hour resolution
omni_data_h = irf_get_data_omni(tint,'n,v,P,Bx,By,Bz,Ms,T,NaNp','omni2');

time = irf_time(omni_data_h(:,1),'epoch>epochtt');
n = irf.ts_scalar(time,omni_data_h(:,2)); n.name = 'n'; n.units = 'cm^{-3}';
v = irf.ts_scalar(time,omni_data_h(:,3)); v.name = 'v'; v.units = 'km s^{-1}';
P = irf.ts_scalar(time,omni_data_h(:,4)); P.name = 'P_{flow}'; P.units = 'nPa';
B = irf.ts_vec_xyz(time,omni_data_h(:,5:7)); B.name = 'B'; B.units = 'nT';
T = irf.ts_scalar(time,omni_data_h(:,9)*K2eV); T.name = 'T'; T.units = 'eV';
nanp = irf.ts_scalar(time,omni_data_h(:,10)); nanp.name = 'n_a/n_p'; 

%Pth = T*n/K2eV*units.eV;
Pth = irf.ts_scalar(time,omni_data_h(:,2)*1e6.*omni_data_h(:,9)*units.kB)*1e9; Pth.name = 'Pth'; Pth.units = 'nPa';

% P - flow pressure
% P (nPa) = (1.67/10**6) * Np*V**2 * (1+ 4*Na/Np)
% for hours with non-fill Na/Np ratios and
% P (nPa) = (2.0/10**6) * Np*V**2
% for hours with fill values for Na/Np 

% min resolution
omni_data_m = irf_get_data_omni(tint,'n,v,P,Bx,By,Bz,Ms,T','omni_min');

time = irf_time(omni_data_m(:,1),'epoch>epochtt');
n_m = irf.ts_scalar(time,omni_data_m(:,2)); n_m.name = 'n'; n_m.units = 'cm^{-3}';
v_m = irf.ts_scalar(time,omni_data_m(:,3)); v_m.name = 'v'; v_m.units = 'km s^{-1}';
P_m = irf.ts_scalar(time,omni_data_m(:,4)); P_m.name = 'P'; P_m.units = 'nPa';
B_m = irf.ts_vec_xyz(time,omni_data_m(:,5:7)); B_m.name = 'B'; B_m.units = 'nT';
T_m = irf.ts_scalar(time,omni_data_m(:,9)*K2eV); T_m.name = 'T'; T_m.units = 'eV';
%nanp_m = irf.ts_scalar(time,omni_data_m(:,10)); nanp_m.name = 'n_a/n_p'; 

%Pth_m = T_m/K2eV*n_m*1e6;
%Pth_m2 = T_m*n_m*1e6*1e9; Pth_m2.name = 'nPa';
Pth_m = irf.ts_scalar(time,omni_data_m(:,2)*1e6.*omni_data_m(:,9)*units.kB)*1e9; Pth_m.name = 'Pth'; Pth_m.units = 'nPa';

%vxyz = irf.ts_scalar(time,omni_data_vxyz(:,2:4)); vxyz.name = 'v'; vxyz.units = 'km s^{-1}';


%%
h = irf_plot(7);

hca = irf_panel('v');
set(hca,'ColorOrder',mms_colors('1x'))
irf_plot(hca,{v,v_m},'comp')
hca.YLabel.String = 'v (km s^{-1})';
hca.YLabel.Interpreter = 'tex';

hca = irf_panel('n');
set(hca,'ColorOrder',mms_colors('1x'))
irf_plot(hca,{n,n_m},'comp')
hca.YLabel.String = 'n (cm^{-3})';
hca.YLabel.Interpreter = 'tex';

if 0 % na/np
  hca = irf_panel('nanp');
  set(hca,'ColorOrder',mms_colors('1x'))
  irf_plot(hca,{nanp},'comp')
  hca.YLabel.String = 'n_a/n_p';
  hca.YLabel.Interpreter = 'tex';
end

hca = irf_panel('T');
set(hca,'ColorOrder',mms_colors('1x'))
irf_plot(hca,{T,T_m},'comp')
hca.YLabel.String = 'T (eV)';
hca.YLabel.Interpreter = 'tex';

hca = irf_panel('Pth');
set(hca,'ColorOrder',mms_colors('1x'))
irf_plot(hca,{Pth,Pth_m},'comp')
hca.YLabel.String = 'P_{therm} (nPa)';

hca = irf_panel('P');
set(hca,'ColorOrder',mms_colors('1x'))
irf_plot(hca,{P,P_m},'comp')
hca.YLabel.String = 'P_{flow} (nPa)';

hca = irf_panel('B');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{B.x,B.y,B.z},'comp')
set(hca,'ColorOrder',mms_colors('xyz'))
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'x','y','z'},[0.98 0.98])
%h = irf_plot({n,v,P,B});


hca = irf_panel('B min');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{B_m.x,B_m.y,B_m.z},'comp')
set(hca,'ColorOrder',mms_colors('xyz'))
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'x','y','z'},[0.98 0.98])
%h = irf_plot({n,v,P,B});

irf_zoom(h,'x',tint)
irf_zoom(h,'y')
