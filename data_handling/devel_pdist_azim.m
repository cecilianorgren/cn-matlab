% Test if the 360 deg spectrogram works
PD_orig = pdist_cleaned;
PD_shift_eye = PD_orig.shift_many([-170 0 0],200,[1 0 0; 0 1 0; 0 0 1],ic);
PD_shift_lmn = PD_orig.shift_many([-170 0 0],200,lmn,ic);
PD_shift_nlm = PD_orig.shift_many([-170 0 0],200,circshift(lmn,1,1),ic);
PD_shift_mnl = PD_orig.shift_many([-170 0 0],200,circshift(lmn,2,1),ic);
PD_shift_mnl_large = PD_orig.shift_many([-2000 0 0],200,circshift(lmn,2,1),ic);
PD_noshift_mnl = PD_orig.shift_many([0 0 0],200,circshift(lmn,2,1),ic);

dir_pitch = irf.ts_vec_xyz(iPDist3.time,repmat(N,iPDist3.length,1));
iPitch_orig  = PD_orig.pitchangles(dir_pitch,12);
iPitch_shift_eye = PD_shift_eye.pitchangles(dir_pitch,12);
iPitch_shift_lmn = PD_shift_lmn.pitchangles(dir_pitch,12);
iPitch_shift_nlm = PD_shift_nlm.pitchangles(dir_pitch,12);
iPitch_shift_mnl = PD_shift_mnl.pitchangles(dir_pitch,12);
iPitch_shift_mnl_large = PD_shift_mnl_large.pitchangles(dir_pitch,12);
iPitch_noshift_mnl = PD_noshift_mnl.pitchangles(dir_pitch,12);

iAzim_orig  = PD_orig.azimuthangle;
iAzim_shift_eye = PD_shift_eye.azimuthangle;
iAzim_shift_lmn = PD_shift_lmn.azimuthangle;
iAzim_shift_nlm = PD_shift_nlm.azimuthangle;
iAzim_shift_mnl = PD_shift_mnl.azimuthangle;
iAzim_shift_mnl_large = PD_shift_mnl_large.azimuthangle;
iAzim_noshift_mnl = PD_noshift_mnl.azimuthangle;

%% Plot
elim = [1500 Inf];

h = irf_plot(8);


if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)    
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1
  hca = irf_panel('azim original');
  specrec = iAzim_orig.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim shifted eye');
  specrec = iAzim_shift_eye.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim shifted lmn');
  specrec = iAzim_shift_lmn.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim shifted nlm');
  specrec = iAzim_shift_nlm.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim shifted mnl');
  specrec = iAzim_shift_mnl.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim shifted mnl large');
  specrec = iAzim_shift_mnl_large.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end
if 1
  hca = irf_panel('azim not shifted mnl');
  specrec = iAzim_noshift_mnl.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.YTick = 0:30:360;
  hca.Layer = 'top';
end

colormap(pic_colors('candy4'))
irf_plot_axis_align(h)

%% Macroparticles

PD = pdist_cleaned;
nt = PD.length;
Vsc = scPot3.resample(PD);
nMP = 5000;
MP = PD.macroparticles('ntot',nMP,'scpot',Vsc);
vL_shift = -170*0;
vL_shift = -999;
E_edges = [PD.ancillary.energy(1,1) - PD.ancillary.delta_energy_minus(1,1), PD.ancillary.energy(1,:) + PD.ancillary.delta_energy_plus(1,:)];
nEnergy = numel(E_edges)-1;
azimuth_edges = -180:10:180;
nAzimuth = numel(azimuth_edges)-1;

dn_tot_all_ = zeros(nt,nAzimuth);
dn_tot_all = zeros(nt,nEnergy,nAzimuth);
df_tot_all = zeros(nt,nEnergy,nAzimuth);
dv_tot_all = zeros(nt,nEnergy,nAzimuth);
for it = 1:nt
  dv = MP(it).dv;
  df = MP(it).df;
  dn = df.*dv;
  vx = MP(it).vx;
  vy = MP(it).vy;
  vz = MP(it).vz;
  vL = [vx vy vz]*L' - vL_shift;
  vM = [vx vy vz]*M';
  vN = [vx vy vz]*N';
  theta_NL = atan2d(vN,vL);
  v2 = (vL.^2 + vM.^2 + vN.^2)*1e6; % (m/s)^2
  E = units.mp*v2/2/units.eV; % eV
  
  [dn_tot_ edges_ mid_ loc_] = histcn([theta_NL],azimuth_edges,'AccumData',dn,'Fun',@sum);
  [dn_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',dn,'Fun',@sum);
  [dv_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',dv,'Fun',@sum);
  [df_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',df,'Fun',@sum);
  % for iloc = 1:nAzimuth    
  %   dv_tot_all(it,iloc) = sum(MP(it).dv(loc==iloc));
  %   df_tot_all(it,iloc) = sum(MP(it).df(loc==iloc));
  % end
  dn_tot_all_(it,:) = dn_tot_;

  dn_tot_all(it,:,:) = dn_tot;
  dv_tot_all(it,:,:) = dv_tot;
  df_tot_all(it,:,:) = df_tot;
  
end

% Not seperated by E
d3v_tot = sum(PD.d3v('mat'),[2:4]);
d3v_tot_az_bin_ = repmat(d3v_tot/nAzimuth,[1,nAzimuth]);
f_tot_2_ = dn_tot_all_./d3v_tot_az_bin_;

% Seperated by E
d3v_tot_per_E = sum(PD.d3v('mat'),[3:4]);
d3v_tot_az_bin = repmat(d3v_tot_per_E/nAzimuth,[1,1,nAzimuth]);
%f_tot_1 = dn_tot_all./dv_tot_all;
f_tot_2 = dn_tot_all./d3v_tot_az_bin;


iAzim_999 = PDist(PD.time,f_tot_2,'azimuthangle',PD.depend{1},mid{2}); % scaling factor to go from 1/m^6 -> 1/cm^6
%iAzim_ = PDist(PD.time,f_tot_1,'azimuthangle',PD.depend{1},mid{2}); % scaling factor to go from 1/m^6 -> 1/cm^6

%% Plot
elim = [1500 Inf];

h = irf_plot(4);


if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)    
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % iAzim_000
  hca = irf_panel('azim iAzim_000');
  specrec = iAzim_000.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.Layer = 'top';
  hca.YTick = -180:30:180;
end
if 1 % iAzim_170
  hca = irf_panel('azim iAzim_170');
  specrec = iAzim_170.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.Layer = 'top';
  hca.YTick = -180:30:180;
end
if 1 % iAzim_999
  hca = irf_panel('azim iAzim_999');
  specrec = iAzim_999.elim(elim).specrec('pitchangle');
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel');  
  hca.Layer = 'top';
  hca.YTick = -180:30:180;
end

colormap('parula')
%colormap(pic_colors('candy4'))
%colormap(flipdim(pic_colors('thermal'),1))
irf_plot_axis_align(h)

c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))