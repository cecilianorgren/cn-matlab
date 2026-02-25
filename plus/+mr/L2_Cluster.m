% mr.L2_Cluster
% Download data
tint = irf.tint('2001-12-03T07:00:00.00Z/2001-12-03T12:00:00.00Z'); %20151112071854
sclist = [1 3 4];
% caa_download(tint,'C?_CP_FGM_FULL',sclist);
% caa_download(tint,'C?_CP_EFW_L2_E3D_INERT',sclist);
% caa_download(tint,'C?_CP_EFW_L2_E3D_GSE',sclist);
% caa_download(tint,'C?_CP_EFW_L2_V3D_INERT',sclist);
% caa_download(tint,'C?_CP_EFW_L2_V3D_GSE',sclist);
% caa_download(tint,'C?_CP_CIS_HIA_ONBOARD_MOMENTS',sclist);
% caa_download(tint,'C?_CP_CIS_CODIF_HS_H1_MOMENTS',sclist);
% caa_download(tint,'C?_CP_PEA_PITCH_SPIN_DEFlux',sclist);

units = irf_units;
% Load data
cd /Users/Cecilia/Data/Cluster/20011203/
c_eval('gseB?=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''ts'');',sclist);
c_eval('dslE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''ts'');',sclist);
c_eval('gseE?=c_caa_var_get(''E_Vec_xyz_GSE__C?_CP_EFW_L2_E3D_GSE'',''ts'');',sclist);
c_eval('dslExB?=c_caa_var_get(''v_drift_ISR2__C?_CP_EFW_L2_V3D_INERT'',''ts'');',sclist);
c_eval('gseExB?=c_caa_var_get(''v_drift_GSE__C?_CP_EFW_L2_V3D_GSE'',''ts'');',sclist);
c_eval('gseVe?=c_caa_var_get(''Data_Velocity_GSE__C?_CP_PEA_MOMENTS'',''ts'');',sclist);
c_eval('gsePe?=c_caa_var_get(''Data_Pressure_GSE__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
c_eval('gsePe? = irf.ts_tensor_xyz(irf_time(gsePe?.t,''epoch>epochtt''),gsePe?.data);',sclist);
c_eval('Ne?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''ts'');',sclist);
c_eval('gseTepar?=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'',''ts'');',sclist);
c_eval('gseTeper?=c_caa_var_get(''Data_Temperature_ComponentPerpendicularToMagField__C?__MOMENTS'',''ts'');',sclist);

c_eval('[caa?,dobj?,pea?,unit?]=c_caa_var_get(''Data__C?_CP_PEA_PITCH_SPIN_DEFlux'');',sclist);

c_eval('gseVi?hia = c_caa_var_get(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''ts'');',sclist);
c_eval('dslVi?hia = c_caa_var_get(''velocity_isr2__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''ts'');',sclist);
c_eval('gsePi?hia = c_caa_var_get(''pressure_tensor__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',sclist);
c_eval('gsePi?hia = irf.ts_tensor_xyz(irf_time(gsePi?hia.t,''epoch>epochtt''),gsePi?hia.data);',[1 3]);
c_eval('Ni?hia = c_caa_var_get(''density__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''ts'');',sclist);
c_eval('gseTipar?hia = c_caa_var_get(''temp_par__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''ts'');',sclist);
c_eval('gseTiper?hia = c_caa_var_get(''temp_perp__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''ts'');',sclist);

c_eval('gseVi?cod = c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''ts'');',sclist);
c_eval('dslVi?cod = c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''ts'');',sclist);
c_eval('gsePi?cod = c_caa_var_get(''pressure__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''mat'');',sclist);
c_eval('gsePi?cod = irf.ts_tensor_xyz(irf_time(gsePi?cod.t,''epoch>epochtt''),gsePi?cod.data);',[1 3]);
c_eval('Ni?cod = c_caa_var_get(''density__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''ts'');',sclist);
c_eval('gseTipar?cod = c_caa_var_get(''T_par__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''ts'');',sclist);
c_eval('gseTiper?cod = c_caa_var_get(''T_perp__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''ts'');',sclist);

% Convective electric field
c_eval('gseVixB? = gseVi?hia.cross(gseB?.resample(gseVi?hia))*1e-3;',sclist)
%% OMNI flux
angles=[7.5000 22.5000 37.5000 52.5000 67.5000 82.5000 97.5000 112.5000 127.5000 142.5000 157.5000 172.5000];
% Summing over pitch angles, with correct weighting due to different solid angle.
c_eval(['eOMNI?.p=squeeze(irf.pitch_angle_average(pea?.data,angles,[7.5 172.5],2));',...
        'eOMNI?.f=pea?.dep_x{2}.data(1,:);',... % energy levels
        'eOMNI?.t=pea?.t;',...
        'eOMNI?.f_label=''Energy [eV]'';',...
        'eOMNI?.p_label={''Electron DEF'',''keV/cm^2-s-str-keV''};'],sclist)
c_eval(['ePA?.f=pea?.dep_x{1}.data;',... % pitch angles
        'ePA?.t=pea?.t;',...
        'ePA?.p=squeeze(nansum(pea?.data,3));',...
        'ePA?.f_label=''Pitchangle [deg]'';',...
        'ePA?.p_label={''Electron DEF'',''keV/cm^2-s-str-keV''};'],sclist)      
                    
return
%% Overview 
ic = 3;

h = irf_plot(5);

hca = irf_panel('gseB');
c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z,gseB?.abs},''comp'')',ic)
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'B_x','B_y','B_z','B_{||}'},[0.98 0.95],'fontsize',14)

if 1 % n from hia
  hca = irf_panel('n');
  c_eval('irf_plot(hca,{gseNi?hia},''comp'')',ic)
  hca.YLabel.String = 'n (cc)';
  irf_legend(hca,{'HIA'},[0.01 0.95],'fontsize',14,'k')
  %irf_legend(hca,{'HIA','COD','PEA'},[0.98 0.95],'fontsize',14)
end

if 0 % n from different instruments
  hca = irf_panel('n');
  c_eval('irf_plot(hca,{gseNi?hia,gseNi?cod,gseNe?},''comp'')',ic)
  hca.YLabel.String = 'n (cc)';
  irf_legend(hca,{'HIA','COD','PEA'},[0.98 0.95],'fontsize',14)
end
if 0 % Ve
  hca = irf_panel('V x');
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'')',ic)
  hca.YLabel.String = 'V_e (km/s)';
end

hca = irf_panel('V hia');
c_eval('irf_plot(hca,{gseVi?hia.x,gseVi?hia.y,gseVi?hia.z},''comp'')',ic)
hca.YLabel.String = 'V_{i} (km/s)';
irf_legend(hca,{'V_x','V_y','V_z'},[0.98 0.95],'fontsize',14)
irf_legend(hca,{'HIA'},[0.01 0.95],'fontsize',14,'k')

hca = irf_panel('gseE');
c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'')',ic)
hca.YLabel.String = 'E (mV/m)';
irf_legend(hca,{'E_x','E_y','E_z'},[0.98 0.95],'fontsize',14)

hca = irf_panel('-gseVixB');
c_eval('irf_plot(hca,{-1*gseVixB?.x,-1*gseVixB?.y,-1*gseVixB?.z},''comp'')',ic)
hca.YLabel.String = '-V_ixB (mV/m)';
irf_legend(hca,{'HIA'},[0.01 0.95],'fontsize',14,'k')
irf_legend(hca,{'E_x','E_y','E_z'},[0.98 0.95],'fontsize',14)

irf_plot_axis_align
irf_zoom(h,'x',tintHT)
%% Plot Rankine-Hugoniot conditions and pressure balance
ic = 3;
units = irf_units;
h = irf_plot(6);

if 0 % gse B
  hca = irf_panel('gseB');
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z,gseB?.abs},''comp'')',ic)
  hca.YLabel.String = 'B (nT)';
  irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.98 0.95],'fontsize',14)
end
if 1 % mva B
  hca = irf_panel('mvaB');
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z,mvaB?.abs},''comp'')',ic)
  hca.YLabel.String = 'B (nT)';
  irf_legend(hca,{'B_L','B_M','B_N','|B|'},[0.98 0.95],'fontsize',14)
end
if 1 % mva Vi hia
  hca = irf_panel('V hia');
  c_eval('irf_plot(hca,{mvaVi?hia.x,mvaVi?hia.y,mvaVi?hia.z},''comp'')',ic)
  hca.YLabel.String = 'V_{i} (km/s)';
  irf_legend(hca,{'V_L','V_M','V_N'},[0.98 0.95],'fontsize',14)
  irf_legend(hca,{'HIA'},[0.01 0.95],'fontsize',14,'k')
end
if 1 % pressure, p and B^2
  hca = irf_panel('pressures');
  c_eval('irf_plot(hca,{gseB?.abs2*1e-18/units.mu0/2*1e9,(gsePe?.resample(gsePi?hia).trace/3+gsePi?hia.trace/3),gseB?.resample(gsePi?hia).abs2*1e-18/units.mu0/2*1e9+(gsePe?.resample(gsePi?hia).trace/3+gsePi?hia.trace/3)},''comp'')',ic)
  hca.YLabel.String = 'Pressure (nPa)';
  irf_legend(hca,{'FGM','PEA+HIA','Total'},[0.98 0.95],'fontsize',14)
end

if 1 % n from hia
  hca = irf_panel('n');
  c_eval('irf_plot(hca,{Ni?hia},''comp'')',ic)
  hca.YLabel.String = 'n (cc)';
  irf_legend(hca,{'HIA'},[0.01 0.95],'fontsize',14,'k')
  %irf_legend(hca,{'HIA','COD','PEA'},[0.98 0.95],'fontsize',14)
end
if 0,   % SUBPLOT: CIS energy spectrogram
     hca = irf_panel('cis')
     irf_plot(hca,'flux__C3_CP_CIS_CODIF_H1_1D_PEF','colorbarlabel',...
         'log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
     %caxis([3.9 7]);
     set(hca,'yscale','log');
     set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
     ylabel(hca,'E_i [eV]');
     irf_legend(hca,{['C' num2str(ic)]},[0.98 0.05],'color','k')
     irf_colormap('default');
 end
 if 1,    % SUBPLOT: PEACE energy spectrogram
     hca = irf_panel('omni e DEFlux')
     irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel',...
         'log10 dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
     caxis([4 9]);
     set(hca,'yscale','log');%,'ylim',[100 3e4]);
     set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])

     irf_legend(hca,{['C' num2str(ic)]},[0.98 0.05],'color','k');
     ylabel('E_e [eV]');
 end
 if 1,    % SUBPLOT: PEACE pitchangle spectrogram
   hca = irf_panel('pa e DEFlux')     
     irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DEFlux','sum_dim2','colorbarlabel',...
         'log10 dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
     caxis([6 9]);     
     set(hca,'ytick',[45 90 135])
     %irf_legend(hca,{['C' num2str(ic)]},[0.98 0.05],'color','k');
     ylabel('Pitchangle [deg.]');
 end
irf_plot_axis_align
irf_zoom(h,'x',tintHT2)
irf_zoom(h(1:4),'y')
hca = 
set(irf_panel('pa e DEFlux'),'ylim',[0 180])

%% Compare velocities, results, use HIA
ic = 1;
npanels = 4;
[h1,h2] = initialize_combined_plot(npanels,3,3,0.4,'vertical');

hca = irf_panel('gseB');
c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z,gseB?.abs},''comp'')',ic)
hca.YLabel.String = 'B (nT)';

hca = irf_panel('Vx');
c_eval('irf_plot(hca,{gseVe?.x,gseVi?hia.x,gseVi?cod.x,gseExB?.resample(gseVi?hia).x},''comp'')',ic)
hca.YLabel.String = 'V_{x} (km/s)';
irf_legend(hca,{'PEA','HIA','COD','EFW'},[0.98 0.95],'fontsize',14)

hca = irf_panel('Vy');
c_eval('irf_plot(hca,{gseVe?.y,gseVi?hia.y,gseVi?cod.y,gseExB?.resample(gseVi?hia).y},''comp'')',ic)
hca.YLabel.String = 'V_{y} (km/s)';
irf_legend(hca,{'PEA','HIA','COD','EFW'},[0.98 0.95],'fontsize',14)

hca = irf_panel('Vz');
c_eval('irf_plot(hca,{gseVe?.z,gseVi?hia.z,gseVi?cod.z,gseExB?.resample(gseVi?hia).z},''comp'')',ic)
hca.YLabel.String = 'V_{z} (km/s)';
irf_legend(hca,{'PEA','HIA','COD','EFW'},[0.98 0.95],'fontsize',14)

%
isub = 1;
% hia efw
hca = h2(isub); isub = isub + 1;
c_eval('scatter(hca,gseVi?hia.x.data,gseExB?.x.resample(gseVi?hia).data)',ic)
hca.XLabel.String = 'v_{HIA}';
hca.YLabel.String = 'v_{EFW}';
irf_legend(hca,'x',[0.1 0.9])

hca = h2(isub); isub = isub + 1;
c_eval('scatter(hca,gseVi?hia.y.data,gseExB?.y.resample(gseVi?hia).data)',ic)
hca.XLabel.String = 'v_{HIA}';
hca.YLabel.String = 'v_{EFW}';
irf_legend(hca,'y',[0.1 0.9])

hca = h2(isub); isub = isub + 1;
c_eval('scatter(hca,gseVi?hia.z.data,gseExB?.z.resample(gseVi?hia).data)',ic)
hca.XLabel.String = 'v_{HIA}';
hca.YLabel.String = 'v_{EFW}';
irf_legend(hca,'z',[0.1 0.9])

% cod efw
hca = h2(isub); isub = isub + 1;
c_eval('scatter(hca,gseVi?cod.x.data,gseExB?.x.resample(gseVi?cod).data)',ic)
hca.XLabel.String = 'v_{COD}';
hca.YLabel.String = 'v_{EFW}';
irf_legend(hca,'x',[0.1 0.9])

hca = h2(isub); isub = isub + 1;
c_eval('scatter(hca,gseVi?cod.y.data,gseExB?.y.resample(gseVi?cod).data)',ic)
hca.XLabel.String = 'v_{COD}';
hca.YLabel.String = 'v_{EFW}';
irf_legend(hca,'y',[0.1 0.9])

hca = h2(isub); isub = isub + 1;
c_eval('scatter(hca,gseVi?cod.z.data,gseExB?.z.resample(gseVi?cod).data)',ic)
hca.XLabel.String = 'v_{COD}';
hca.YLabel.String = 'v_{EFW}';
irf_legend(hca,'z',[0.1 0.9])

% hia cod
hca = h2(isub); isub = isub + 1;
c_eval('scatter(hca,gseVi?hia.x.resample(gseVi?cod).data,gseVi?cod.x.data)',ic)
hca.XLabel.String = 'v_{HIA} (km/s)';
hca.YLabel.String = 'v_{COD} (km/s)';
irf_legend(hca,'x',[0.1 0.9])

hca = h2(isub); isub = isub + 1;
c_eval('scatter(hca,gseVi?hia.y.resample(gseVi?cod).data,gseVi?cod.y.data)',ic)
hca.XLabel.String = 'v_{HIA} (km/s)';
hca.YLabel.String = 'v_{COD} (km/s)';
irf_legend(hca,'y',[0.1 0.9])

hca = h2(isub); isub = isub + 1;
c_eval('scatter(hca,gseVi?hia.z.resample(gseVi?cod).data,gseVi?cod.z.data)',ic)
hca.XLabel.String = 'v_{HIA} (km/s)';
hca.YLabel.String = 'v_{COD} (km/s)';
irf_legend(hca,'z',[0.1 0.9])

xlim = [-300 300];
ylim = [-300 300];
for ip = 1:9;
  hold(h2(ip),'on')
  plot(h2(ip),xlim,ylim)
  hold(h2(ip),'off')
  h2(ip).XLim = xlim;
  h2(ip).YLim = ylim;
end
%% Do Walen analysis
tintHT0 = irf.tint('2001-12-03T10:57:40.000Z/2001-12-03T10:59:00.000Z'); 
tintHT1 = irf.tint('2001-12-03T11:20:26.953Z/2001-12-03T11:25:30.804Z'); % later
tintHT1 = irf.tint('2001-12-03T11:23:00.000Z/2001-12-03T11:25:30.804Z'); % later
%tintHT1 = irf.tint('2001-12-03T11:23:00.000Z/2001-12-03T11:30:30.804Z'); % later wider
tintHT2_ = irf.tint('2001-12-03T10:57:25.000Z/2001-12-03T10:59:30.000Z'); % Retino outbound
tintHT2 = irf.tint('2001-12-03T10:57:29.266Z/2001-12-03T10:59:27.177Z'); % Retino outbound
tintHT2 = irf.tint('2001-12-03T10:57:53.000Z/2001-12-03T10:59:27.177Z'); % Retino outbound
tintWa2  = irf.tint('2001-12-03T10:57:55.000Z/2001-12-03T10:58:40.000Z'); % Retino outbound
tintHT = tintHT1;
tintHT_ = tintHT2_;
tintWa = tintWa2;
[mvaB3,l,lmn]=irf_minvar(gseB3.tlim(tintHT));
mvaVi3hia = gseVi3hia*lmn';
mvaB3 = gseB3*lmn';
%%
% Find HT frame
[vht,eht,dvht,p,cc] = irf_vht(-gseVixB3.tlim(tintHT)*lmn',gseB3.tlim(tintHT)*lmn');
cc_ht = cc(1,2);

% Do Walen test, scatter plot
ic = 3;
%[Fpe,Fce,Fuh,Fpp,Fcp,FpO,FcO,Va,Vte,Le] = irf_plasma_calc(gseB3,gseNi3hia,gseNi3hia*0,gseTepar3,gseTipar4hia,1);
c_eval('VA?     = 1e-9*mvaB?.resample(Ni?hia).abs/sqrt(units.mu0*(Ni?hia*1e6)*units.mp)*1e-3;',ic);
c_eval('mvaVA?  = 1e-9*mvaB?.resample(Ni?hia)/sqrt(units.mu0*(Ni?hia*1e6)*units.mp)*1e-3;',ic);
c_eval('v = mvaVi?hia;',ic)
c_eval('va = mvaVA?.resample(v);',ic)
vvht = v-1*vht;

npanels = 3;
[h1,h2] = initialize_combined_plot(npanels,2,1,0.6,'vertical');

isub = 1;
hca = h1(isub); isub = isub + 1;
irf_plot(hca,mvaB3)
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'L','M','N'},[0.98 0.95],'fontsize',14)
  
hca = h1(isub); isub = isub + 1;
irf_plot(hca,vvht)
ylabel(hca,'v-v_{HT} (km/s)','interpreter','tex')
irf_legend(hca,{'L','M','N'},[0.98 0.95],'fontsize',14)
  
hca = h1(isub); isub = isub + 1;
irf_plot(hca,va)
ylabel(hca,'v_{A} (km/s)','interpreter','tex')
irf_legend(hca,{'L','M','N'},[0.98 0.95],'fontsize',14)
irf_legend(hca,{'projected on B'},[0.98 0.1],'fontsize',16,'color',[0 0 0])


for ip = 1:3;
  h1(ip).FontSize = 16;
end
for ip = 2:3;
  h1(ip).YLim = [-400 400];
end
irf_zoom(h1(1),'y')
hmarkWa = irf_pl_mark(h1,tintWa,'y'); for im = 1:numel(hmarkWa), hmarkWa(im).FaceAlpha = 1; end
hmarkHT = irf_pl_mark(h1,tintHT,'g'); for im = 1:numel(hmarkHT), hmarkHT(im).FaceAlpha = 0.5; end
irf_plot_axis_align
irf_zoom(h1,'x',tintHT_)


h2(1).Position(2) = h2(1).Position(2)-0.05;
h2(1).Position(1) = h2(1).Position(1)-0.05;
h2(2).Position(1) = h2(2).Position(1)-0.05;
%%
isub = 1;
colors = mms_colors('matlab');
hca = h2(isub); isub = isub + 1;
% info string on HT analysis
strvht=['V_{HT}=' num2str(irf_abs(vht,1),3) ' [ ' num2str(irf_norm(vht),' %5.2f') '] km/s GSE'];
strdvht=['\Delta V_{HT}=' num2str(irf_abs(dvht,1),3) ' [ ' num2str(irf_norm(dvht),' %5.2f') '] km/s GSE'];
strint = [irf_time(tintHT(1).epochUnix,'utc_yyyy-mm-dd HH:MM:SS.mmm') ' - ' irf_time(tintHT(2).epochUnix,'utc_yyyy-mm-dd HH:MM:SS.mmm')];
strsl  = ['slope=' num2str(p(1),3) '  offs=' num2str(p(2),2)];
stroffs= ['cc=' num2str(cc(1,2),3)];
%irf_vht_plot(-gseVixB3.tlim(tintHT)*lmn',gseB3.tlim(tintHT)*lmn');

e = -gseVixB3.tlim(tintHT)*lmn'; e = e.data;
%plot(hca,eht.data(:,1),e(:,1),'o',eht.data(:,2),e(:,2),'o',eht.data(:,3),e(:,3),'o');
plot(hca,eht.data(:,1),e(:,1),'o',eht.data(:,2),e(:,2),'o');
axes(hca)

legend(hca,'L','M','N','location','southeast')
hca.FontSize = 16;
hold(hca,'on')
hca.XLim = 4*[-1 1];
hca.YLim = 4*[-1 1];
xlim = hca.XLim;
ylim = hca.YLim;
plot(hca,1e3*[-1 1],1e3*[-1 1],'k-')
plot(hca,1e3*[-1 1],p(1)*1e3*[-1 1]+p(2),'--','color',0.6*[1 1 1])
%text(-350,100,sprintf(' slope = %.2f \n offset = %.0f km/s \n cc = %.2f',fitx(1),fitx(2),ccx),'fontsize',16,'color',colors(1,:))
hca.XLim = xlim;
hca.YLim = ylim;
hold(hca,'off')
axis(hca,'equal')
axis(hca,'square')
hca.Box = 'on';
hca.XLabel.String = 'E (mV/m)';
hca.YLabel.String = '-V_ixB (mV/m)';
text(hca.XLim(1),hca.YLim(2),9,sprintf(' %s \n %s \n %s \n %s \n %s',strint,strvht,strdvht,strsl,stroffs),'fontsize',12,'color',[0 0 0],'verticalalignment','bottom','horizontalalignment','left')

hca = h2(isub); isub = isub + 1;
ccx = xcorr(vvht.tlim(tintWa).x.data,va.tlim(tintWa).x.data,0,'coeff');
ccy = xcorr(vvht.tlim(tintWa).y.data,va.tlim(tintWa).y.data,0,'coeff');
ccz = xcorr(vvht.tlim(tintWa).z.data,va.tlim(tintWa).z.data,0,'coeff');

fitx = polyfit(vvht.tlim(tintWa).x.data,va.tlim(tintWa).x.data,1);
fity = polyfit(vvht.tlim(tintWa).y.data,va.tlim(tintWa).y.data,1);
fitz = polyfit(vvht.tlim(tintWa).z.data,va.tlim(tintWa).z.data,1);

scatter(hca,vvht.tlim(tintWa).x.data,va.tlim(tintWa).x.data); hold(hca,'on')
scatter(hca,vvht.tlim(tintWa).y.data,va.tlim(tintWa).y.data); 
scatter(hca,vvht.tlim(tintWa).z.data,va.tlim(tintWa).z.data); 
hca.XLim = 400*[-1 1];
hca.YLim = 400*[-1 1];
xlim = hca.XLim;
ylim = hca.YLim;
plot(hca,1e3*[-1 1],1e3*[-1 1],'k-')
plot(hca,1e3*[-1 1],fitx(1)*1e3*[-1 1]+fitx(2),'--','color',colors(1,:))
axes(hca)
text(hca.XLim(1),hca.YLim(2),sprintf(' slope = %.2f \n offset = %.0f km/s \n cc = %.2f',fitx(1),fitx(2),ccx),'fontsize',16,'color',colors(1,:),'verticalalignment','top','horizontalalignment','left')
hca.XLim = xlim;
hca.YLim = ylim;
hold(hca,'off')
legend(hca,'L','M','N','location','southeast')

hca.XLabel.String = 'v-v_{HT} (km/s)';
hca.YLabel.String = 'v_{A} (km/s)';
axis(hca,'equal')
axis(hca,'square')
hca.Box = 'on';

hca.FontSize = 16;


