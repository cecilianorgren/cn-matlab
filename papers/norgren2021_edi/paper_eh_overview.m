%% Load data
%ic = 1:4;
units = irf_units;
tint_brst = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint = tint_brst + [+5 -5]; % using the above edges causes problem with new EDI files because they have different versions that adjoining file
tint_figure = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.630Z');
tint = tint_figure + [-5 5];
localuser = datastore('local','user');
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
db_info = datastore('mms_db');   

c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',1:4);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); ',1:4);

c_eval('gseR? = mms.get_data(''R_gse'',tint_brst,?); gseR? = gseR?.resample(tint_figure(1))',1:4)

c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',1:4);

%c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
%c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);

c_eval('ePitch?_flux_edi = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',1:4)

c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+1));',1:4)
c_eval('ePDist? = ePDist?.tlim(tint);',1:4)

%% Prepare data
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',1:4)

n0 = 0.04; % for normalization of density perturbation later

vph = -8500e3;
c_eval('[phi?,phi_progressive?,phi_ancillary?] = get_phi(gseE?par,vph,tint_zoom,tint_zoom);',1:4)

c_eval('vtrap? = irf.ts_scalar(phi?.time,sqrt(units.e*phi?.data/units.me)*1e-3); vtrap?.units = ''km/s'';',1:4)
%c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = ''km/s'';',1:4) % km/s
c_eval('gseVExB? = cross(gseE?,gseB?.resample(gseE?.time))/gseB?.resample(gseE?.time).abs/gseB?.resample(gseE?.time).abs*1e3; gseVExB?.units = ''km/s'';',1:4) % km/s
c_eval('wExB? = irf.ts_scalar(gseVExB?.time,units.me*gseVExB?.abs2.data*10^6/2/units.eV); wExB?.units = ''eV'';',1:4) % eV

% Density pertubratin from 4 sc
c_eval('R? = gseR?.resample(gseE1)*1e3;',1:4)
c_eval('E? = gseE?.resample(gseE1)*1e-3;',1:4)
[divE,avE]=c_4_grad('R?','E?','div'); divE.name = 'div E';
dn_divE = divE*units.eps0/units.e*1e-6; 
dn_divE.name = 'dn from div E';
dn_divE.units = 'cm^-3';
  
% Diff E, proxy of Gauss law density perturbation
c_eval('dtE? = gseE?.time(2) - gseE?.time(1);',1:4)
c_eval('diffE? = irf.ts_vec_xyz([gseE?.time(1:end-1)+0.5*dtE?],diff(gseE?.data,1));',1:4)
c_eval('diffE?par = diffE?.dot(gseB?.resample(diffE?).norm);',1:4)
c_eval('dn_E?par = diffE?par*1e-3/dtE?*units.eps0/units.e;',1:4)

% Prepare flux in closest FPI energy channel.
if 1 % Make the data stepfunction-like, so that it shows the accumulation time
  %%
  c_eval('ePitch?_flux_edi_apar = ePitch?_flux_edi.palim(180).resample(ePDist?);',1:4);
  
  c_eval('ePitch?_fpi = ePDist?.pitchangles(dmpaB?,16);',1:4); 
  c_eval('ePitch?_flux_fpi = ePitch?_fpi.flux;',1:4);
  c_eval('ePitch?_flux_fpi_apar500 = ePitch?_flux_fpi.elim(500).palim(180);',1:4);
  c_eval('ePitch?_flux_fpi_par500 = ePitch?_flux_fpi.elim(500).palim(0);',1:4);
  
  c_eval('ePitch?_fpi2 = ePDist?.pitchangles(dmpaB?,180-[22.5 0]);',1:4);
  c_eval('ePitch?_flux_fpi2 = ePitch?_fpi2.flux;',1:4);
  c_eval('ePitch?_flux_fpi2_apar500 = ePitch?_flux_fpi2.elim(500).palim(180);',1:4);  
  % Make the data stepfunction-like, so that it shows the accumulation time
  dt = 0.015;
  for iic = 1:4
    c_eval('pitch_tmp = ePitch?_flux_edi_apar;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_edi_apar_step = irf.ts_scalar(newtime,newdata);',iic)
    
    
    c_eval('pitch_tmp = ePitch?_flux_fpi2_apar500;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_fpi2_apar500_step = irf.ts_scalar(newtime,newdata);',iic)
    
    c_eval('pitch_tmp = ePitch?_flux_fpi_apar500;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_fpi_apar500_step = irf.ts_scalar(newtime,newdata);',iic)
    
    c_eval('pitch_tmp = ePitch?_flux_fpi_par500;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_fpi_par500_step = irf.ts_scalar(newtime,newdata);',iic)
  end
end

%% Prepare data, potential
% First run
% mms_20170706_135303.load_data
% mms_20170706_135303.prepare_data_single_sc
% mms_20170706_135303.prepare_data_multi_sc
tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); % if showing 4 sc epar
ne1_mean = mean(ne1.tlim(tint_zoom).data);
n0 = 0.04;
if 0 % dont need reduced distribution for this plot
  eDist = ePDist1.tlim(tint);
  % remove background
  nSecondary = [5];
  nPhoto = 0;
  %[eDist_nobg] = mms.remove_edist_background(eDist_orig);
  c_eval('[eDist_nobg?] = mms.remove_edist_background(eDist,''nSecondary'',nSecondary(?),''Nphotoe_art'',nPhoto,''ZeroNaN'',0);',1:numel(nSecondary))

  eint = [000 40000];
  vint = [-Inf Inf];
  scpot = scPot1.resample(eDist);
  lowerelim = scpot*0 + 00;
  %tic; ef1D = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'vint',vint,'scpot',scpot,'lowerelim',lowerelim); toc % reduced distribution along B
  tic; ef1D = eDist_nobg1.reduce('1D',dmpaB1.resample(eDist).norm,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'nMC',500); toc % reduced distribution along B
end
vph = -8500e3;
c_eval('[phi?,phi_progressive?,phi_ancillary?] = get_phi(gseE?par,vph,tint_zoom,tint_zoom);',1:4)

% Prepare flux in closest FPI energy channel.
if 1 % Make the data stepfunction-like, so that it shows the accumulation time
  %%
  c_eval('ePitch?_flux_edi_apar = ePitch?_flux_edi.palim(180).resample(ePDist?);',1:4);
  
  c_eval('ePitch?_fpi = ePDist?.pitchangles(dmpaB?,16);',1:4); 
  c_eval('ePitch?_flux_fpi = ePitch?_fpi.flux;',1:4);
  c_eval('ePitch?_flux_fpi_apar500 = ePitch?_flux_fpi.elim(500).palim(180);',1:4);
  c_eval('ePitch?_flux_fpi_par500 = ePitch?_flux_fpi.elim(500).palim(0);',1:4);
  
  c_eval('ePitch?_fpi2 = ePDist?.pitchangles(dmpaB?,180-[22.5 0]);',1:4);
  c_eval('ePitch?_flux_fpi2 = ePitch?_fpi2.flux;',1:4);
  c_eval('ePitch?_flux_fpi2_apar500 = ePitch?_flux_fpi2.elim(500).palim(180);',1:4);  
  % Make the data stepfunction-like, so that it shows the accumulation time
  dt = 0.015;
  for iic = 1:4        
    c_eval('pitch_tmp = ePitch?_flux_edi_apar;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_edi_apar_step = irf.ts_scalar(newtime,newdata);',iic)
    
    
    c_eval('pitch_tmp = ePitch?_flux_fpi2_apar500;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_fpi2_apar500_step = irf.ts_scalar(newtime,newdata);',iic)
    
    c_eval('pitch_tmp = ePitch?_flux_fpi_apar500;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_fpi_apar500_step = irf.ts_scalar(newtime,newdata);',iic)
    
    c_eval('pitch_tmp = ePitch?_flux_fpi_par500;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_fpi_par500_step = irf.ts_scalar(newtime,newdata);',iic)
  end
end

%% Plot, zoom only
ic = 1;
npanels = 8;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
fontsize = 12;

v_for_density_scaling = 8500e3; % m/s

tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z')+0*[-8 8];
tint_zoom = irf.tint('2017-07-06T13:54:05.30Z/2017-07-06T13:54:05.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.40Z/2017-07-06T13:54:06.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); % if showing 4 sc epar
tint_zoom = phi1.time([1 end]);

% load eh data
manual = edi_event_manual_dt;

% data_tmp = load(sprintf('/Users/%s/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat',localuser));
% obs_eh_properties = data_tmp.EH_properties;
% obs_lpp = obs_eh_properties.Lpp; % peak to peak length
% obs_potential = obs_eh_properties.potential;
% obs_vtrap = sqrt(2*units.e*obs_potential/units.me)*1e-3; % km/s
% obs_potential_max = obs_eh_properties.max_potential;
% obs_velocity = obs_eh_properties.vel;
% obs_neh = numel(obs_velocity);
% c_eval('obs_t0_epoch_mms? = obs_eh_properties.time_mms?;')
% c_eval('obs_phi? = irf.ts_scalar(obs_t0_epoch_mms?,obs_potential(:,?));')
% c_eval('obs_vph? = irf.ts_scalar(obs_t0_epoch_mms?,obs_velocity);')
% c_eval('obs_vtrap? = irf.ts_scalar(obs_t0_epoch_mms?,obs_vtrap(:,?));')
% c_eval('obs_vph?.data(isnan(obs_potential(:,?))) = NaN;')

% common time shifts
dt = [0.0000  -0.0012  -0.0009  -0.0012]; % old
dt = [0   -0.0010   -0.0010   -0.0013]; % new dt = mean(cat(1,manual.dt),1);
dt0 = 0.0008;
dt = dt + dt0;

flow = 3;

if 1 % E par, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'MMS 1';'MMS 2';'MMS 3';'MMS 4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % E perp, abs, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp abs dt');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseE1perp.abs,gseE2perp.abs,gseE3perp.abs,gseE4perp.abs},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  
  hca.YLabel.String = {'|E_{\perp}|','(mV/m)'};  
end
if 1 % E perp, abs, filt, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp abs dt filt');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseE1perp.filt(flow,0,[],5).abs,gseE2perp.filt(flow,0,[],5).abs,gseE3perp.filt(flow,0,[],5).abs,gseE4perp.filt(flow,0,[],5).abs},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  format_ms = '%.1f';
  irf_legend(hca,{sprintf('f > %g Hz',flow)},[0.05 0.98],'fontsize',12);
  hca.YLabel.String = {'|E_{\perp}|','(mV/m)'};  
end
if 1 % Phi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{phi1,phi2,phi3,phi4},'comp','dt',dt);  
  hca.YLabel.String = {'\phi','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end
if 0 % ExB energy
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('wExB');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{wExB1.filt(flow,0,[],5),wExB2.filt(flow,0,[],5),wExB3.filt(flow,0,[],5),wExB4.filt(flow,0,[],5)},'comp','dt',dt);  
  hca.YLabel.String = {'W_{ExB}','(eV)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')  
end
if 0 % vtrap
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('vtrap');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{vtrap1,vtrap2,vtrap3,vtrap4},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))    
  hca.YLabel.String = {'v_{trap}','(km/s)'};  
end
if 0 % |vExB|
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('v ExB abs');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseVExB1.filt(flow,0,[],5).abs,gseVExB2.filt(flow,0,[],5).abs,gseVExB3.filt(flow,0,[],5).abs,gseVExB4.filt(flow,0,[],5).abs},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))    
  hca.YLabel.String = {'|v_{ExB}|','(km/s)'};  
  format_ms = '%.1f';
  irf_legend(hca,{sprintf('f > %g Hz',flow)},[0.05 0.98],'fontsize',12);  
end
if 1 % charge perturbation
  sub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('density perturbation');
  set(hca,'ColorOrder',mms_colors('1234b'))
  % cm scale 
  units_scaling = 1e-3;
  hh = irf_plot(hca,{dn_E1par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E2par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E3par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E4par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_divE/units_scaling,...
                     },'comp','dt',[dt 0]); 
  hh(1).Children(1).LineWidth = 2;
  hca.YLabel.String = {'\delta n',sprintf('(10^{%g} cm^{-3})',log10(units_scaling))};
  hca.YLabel.FontSize = fontsize;
  ylabel(hca,hca.YLabel.String,'interpreter','tex')

  if 1 % 4sc
    set(hca,'ColorOrder',mms_colors('b'))
    hleg_4sc = irf_legend(hca,{'4 sc'},[0.99 0.98],'fontsize',12);
    hleg_4sc.FontWeight = 'bold';
  elseif 1 % legend mms1 mms2 mms43 mms4, 4sc
    set(hca,'ColorOrder',mms_colors('1234b'))
    irf_legend(hca,{'mms1';'mms2';'mms3';'mms4';'4 sc'},[1.02 0.9],'fontsize',12);
  else
    set(hca,'ColorOrder',mms_colors('b'))
    irf_legend(hca,{'4 sc'},[1.02 0.9],'fontsize',12);
  end
end
if 1 % charge perturbation
  sub = isub + 1;
  zoomy = [zoomy isub];  
  hca = irf_panel('scpot filt');
  %ffilt = 100;
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{scPot1.filt(ffilt,0,[],5),scPot2.filt(ffilt,0,[],5),scPot3.filt(ffilt,0,[],5),scPot4.filt(ffilt,0,[],5)},'comp','dt',dt)
  hca.YLabel.String = {'\delta V_{sc}','(V)'};
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,sprintf('f>%g Hz',ffilt),[0.05 0.98])
end
if 1 % edi flux 180 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  palim = [168 180];  
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6,ePitch3_flux_edi.palim(palim)*1e-6,ePitch4_flux_edi.palim(palim)*1e-6},'comp','dt',dt);
  %irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6*0.4,ePitch3_flux_edi.palim(palim)*1e-6*0.5,ePitch4_flux_edi.palim(palim)*1e-6*0.5},'comp','dt',dt);  
  hca.YLabel.String = {'j_e^{EDI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{EDI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'\theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
end
if 1 % edi flux 180 4sc, averaged to fpi and stepped
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux apar step');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_edi_apar_step*1e-6,ePitch2_flux_edi_apar_step*1e-6,ePitch3_flux_edi_apar_step*1e-6,ePitch4_flux_edi_apar_step*1e-6},'comp')  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{EDI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
  irf_legend(hca,{'\theta = [168.25, 180]^o'},[0.05 0.99],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'no time shift'},[0.05 0.8],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'averaged to FPI timeline'},[0.98 0.98],'fontsize',12,'color',[0 0 0]);
end
if 1 % fpi flux 180 4sc, 22.50
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi flux apar 2');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_fpi2_apar500_step*1e-6,ePitch2_flux_fpi2_apar500_step*1e-6,ePitch3_flux_fpi2_apar500_step*1e-6,ePitch4_flux_fpi2_apar500_step*1e-6},'comp')  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{FPI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
  irf_legend(hca,{'\theta = [157.50, 180]^o'},[0.05 0.99],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'no time shift'},[0.05 0.8],'fontsize',12,'color',[0 0 0]);
end


%c_eval('h(?).Position(2) = h(?).Position(2)-0.03;',isub_short)
%hcb_fred.Position(2) = hcb_fred.Position(2)+0.02;
%c_eval('h(?).Position(2) = h(?).Position(2)+0.02;',isub_long)


irf_zoom(h,'x',phi1.time)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends_color = {'k','k','k','k','k','k','k','k','k','k','k','k'};
for ipanel = 1:npanels
  irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  h(ipanel).YLabel.FontSize = fontsize;
  h(ipanel).FontSize = fontsize;
end

hca = irf_panel('edi flux'); hca.YLim = [0 7.999];
hca = irf_panel('edi flux apar step'); hca.YLim = [0 7.999];
hca = irf_panel('fpi flux apar 2'); hca.YLim = [0 7.999];
hca = irf_panel('E par dt'); hca.YLim = [-70 60];

% Give the two edi panels the same ylabel
hca = irf_panel('edi flux'); pos1 = hca.YLabel.Position;
hca.YLabel.Position(2) = 0;
hca = irf_panel('edi flux apar step'); pos2 = hca.YLabel.Position;
hca.YLabel = [];

%c_ev

doDoubleAxis = 1; % dn
if doDoubleAxis  
  hca = irf_panel('density perturbation');
  ax1 = hca;
  ax2 = axes('Position',get(ax1,'Position'));
  ax2.YLim = ax1.YLim*1e-3/n0;    
  set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required
  set(ax2,'YAxisLocation','right');
  set(ax2,'Color','none','box','off'); % color of axis      
  %set(ax2,'XColor','k','YColor',colors(2,:)); % color of axis lines and numbers
  %set(ax1,'XColor','k','YColor',colors(1,:)); % color of axis lines and numbers
  irf_timeaxis(ax2,'nolabels')
  ax2.XLabel.String = [];
  ax2.YLabel.String = {'\delta n/n'};
  ax2.YLabel.Interpreter = 'tex';    
  ax2.YTick = hca.YTick*1e-3/n0;  
  ax2.FontSize = fontsize;
end  
  
%irf_plot_axis_align(h)
%ax.YLabel.Position(1) = 1.07;

set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required
