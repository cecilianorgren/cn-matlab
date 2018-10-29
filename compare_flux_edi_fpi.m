%% Load data
ic = 1:4;
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint = irf.tint('2015-11-12T07:18:54.00Z/2015-11-12T07:19:45.00Z');
localuser = datastore('local','user');
mms.db_init('local_file_db','/Volumes/Nexus/data');
db_info = datastore('mms_db');   

% Magnetic field
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);

% Spacecraft potential
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);

% Electron distributions
%c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('tic; [ePDist?,~] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)

% Electron density
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);

%% Flux EDI
if 0
  node = 1;
  c_eval('flux0_mms? = mms.db_get_ts(''mms?_edi_brst_l2_amb-pm2'',''mms?_edi_flux!_0_brst_l2'',tint);',ic,node)
  c_eval('flux180_mms? = mms.db_get_ts(''mms?_edi_brst_l2_amb-pm2'',''mms?_edi_flux!_180_brst_l2'',tint);',ic,node)
else
  node = 2:3;
  c_eval('flux0_node!_mms? = mms.db_get_ts(''mms?_edi_brst_l2_amb'',''mms?_edi_flux!_0_brst_l2'',tint);',ic,node)
  c_eval('flux0_mms? = [];',ic)
  c_eval('flux0_mms? = [flux0_mms? flux0_node!_mms?.data];',ic,node)
  c_eval('flux0_mms? = irf.ts_scalar(flux0_node!_mms?.time,mean(flux0_mms?,2));',ic,node(1))
  
  c_eval('flux180_node!_mms? = mms.db_get_ts(''mms?_edi_brst_l2_amb'',''mms?_edi_flux!_180_brst_l2'',tint);',ic,node)
  c_eval('flux180_mms? = [];',ic)
  c_eval('flux180_mms? = [flux180_mms? flux180_node!_mms?.data];',ic,node)
  c_eval('flux180_mms? = irf.ts_scalar(flux180_node!_mms?.time,mean(flux180_mms?,2));',ic,node(1))
end
%%
% Prepare FPI flux
%edi_pitchangles = [0:11.25:180];
edi_pitchangles = [180-11.25 180];

c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,edi_pitchangles);',1:4)
c_eval('eFlux? = ePitch?.flux;',1:4)

%%
units = irf_units;

% EDI energy and corresponding velocity
E_edi = 250; % eV
v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
dE_edi = 0.05*E_edi; % eV

E_edi_plus = E_edi + dE_edi;
E_edi_minus = E_edi - dE_edi;
v_edi_plus = sqrt(2*units.e*E_edi_plus./units.me); % m/s
v_edi_minus = sqrt(2*units.e*E_edi_minus./units.me); % m/s
dv_edi_minus = v_edi_minus - v_edi;
dv_edi_plus = v_edi_plus - v_edi;
dv_edi = dv_edi_plus - dv_edi_minus; % m/s

% FPI
E_fpi = E_edi;
v_fpi = sqrt(2*units.e*E_fpi/units.me); % m/s
ePDist_elim = ePDist1.elim(E_fpi);

E_fpi_plus = E_fpi + ePDist_elim.ancillary.delta_energy_plus(1,1);
E_fpi_minus = E_fpi - ePDist_elim.ancillary.delta_energy_minus(1,1);
v_fpi_plus = sqrt(2*units.e*E_fpi_plus./units.me); % m/s
v_fpi_minus = sqrt(2*units.e*E_fpi_minus./units.me); % m/s
dv_fpi_minus = v_fpi_minus - v_fpi;
dv_fpi_plus = v_fpi_plus - v_fpi;
dv_fpi = dv_fpi_plus - dv_fpi_minus; % m/s




dv_edi/dv_fpi
%% Figure, tseries of edi and fpi flux
npanels = 7;
%h = irf_plot(npanels);
[h, h2] = initialize_combined_plot(npanels,4,1,0.7,'vertical');

isub = 1;
flux_scale = 1e-6;

if 1 % Flux EDI
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{flux180_mms1*flux_scale,flux180_mms2*flux_scale,flux180_mms3*flux_scale,flux180_mms4*flux_scale},'comp')
  hca.YLabel.String = {'Flux',sprintf('(10^%g cm^{-1}s^{-1}sr^{-1})',log10(1/flux_scale))};
  irf_legend(hca,'EDI',[0.02 0.98],[0 0 0])
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.01 0.98],[0 0 0])
end
if 1 % Flux EDI
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{flux180_mms1.resample(ePitch1)*flux_scale,flux180_mms2.resample(ePitch2)*flux_scale,flux180_mms3.resample(ePitch3)*flux_scale,flux180_mms4.resample(ePitch4)*flux_scale},'comp')
  hca.YLabel.String = {'Flux',sprintf('(10^%g cm^{-1}s^{-1}sr^{-1})',log10(1/flux_scale))};
  irf_legend(hca,'EDI resampled to FPI',[0.02 0.98],[0 0 0])
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.01 0.98],[0 0 0])
end
if 1 % Flux FPI
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{ePitch1.elim(500).flux*flux_scale,ePitch2.elim(500).flux*flux_scale,ePitch3.elim(500).flux*flux_scale,ePitch4.elim(500).flux*flux_scale},'comp')
  hca.YLabel.String = {'Flux',sprintf('(10^%g cm^{-1}s^{-1}sr^{-1})',log10(1/flux_scale))};
  irf_legend(hca,'FPI',[0.02 0.98],[0 0 0])
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.01 0.98],[0 0 0])
end
if 1
  ic = 1;
  hca = h(isub); isub = isub + 1;
  c_eval('irf_plot(hca,{flux180_mms?*flux_scale,ePitch?.elim(500).flux*dv_edi/dv_fpi*flux_scale},''comp'')',ic)
  hca.YLabel.String = {'Flux',sprintf('(10^%g cm^{-1}s^{-1}sr^{-1})',log10(1/flux_scale))};
  irf_legend(hca,{'EDI';'FPI'},[1.01 0.98],[0 0 0])
  irf_legend(hca,sprintf('mms%g',ic),[0.02 0.98],[0 0 0])
end
if 1
  ic = 2;
  hca = h(isub); isub = isub + 1;
  c_eval('irf_plot(hca,{flux180_mms?*flux_scale,ePitch?.elim(500).flux*dv_edi/dv_fpi*flux_scale},''comp'')',ic)
  hca.YLabel.String = {'Flux',sprintf('(10^%g cm^{-1}s^{-1}sr^{-1})',log10(1/flux_scale))};
  irf_legend(hca,{'EDI';'FPI'},[1.01 0.98],[0 0 0])
  irf_legend(hca,sprintf('mms%g',ic),[0.02 0.98],[0 0 0])
end
if 1
  ic = 3;
  hca = h(isub); isub = isub + 1;
  c_eval('irf_plot(hca,{flux180_mms?*flux_scale,ePitch?.elim(500).flux*dv_edi/dv_fpi*flux_scale},''comp'')',ic)
  hca.YLabel.String = {'Flux',sprintf('(10^%g cm^{-1}s^{-1}sr^{-1})',log10(1/flux_scale))};
  irf_legend(hca,{'EDI';'FPI'},[1.01 0.98],[0 0 0])
  irf_legend(hca,sprintf('mms%g',ic),[0.02 0.98],[0 0 0])
end
if 1
  ic = 4;
  hca = h(isub); isub = isub + 1;
  c_eval('irf_plot(hca,{flux180_mms?*flux_scale,ePitch?.elim(500).flux*dv_edi/dv_fpi*flux_scale},''comp'')',ic)
  hca.YLabel.String = {'Flux',sprintf('(10^%g cm^{-1}s^{-1}sr^{-1})',log10(1/flux_scale))};
  irf_legend(hca,{'EDI';'FPI'},[1.01 0.98],[0 0 0])
  irf_legend(hca,sprintf('mms%g',ic),[0.02 0.98],[0 0 0])
end

h(1).Title.String = 'Comparison of flux for EDI and FPI at 180^\circ';
irf_plot_axis_align
irf_zoom(h,'x',ePitch1.time)

for ipanel = 1:npanels
  h(ipanel).FontSize = 10;
end

isub = 1;
for ic = 1:4
  if 1
    hca = h2(isub); isub = isub + 1;    
    c_eval('scatter(hca,ePitch?.elim(500).flux.data*dv_edi/dv_fpi*flux_scale,flux180_mms?.resample(ePitch?).data*flux_scale,''.'')',ic)
    hca.XLabel.String = sprintf(' Flux FPI (10^%g cm^{-1}s^{-1}sr^{-1})',log10(1/flux_scale));
    hca.YLabel.String = sprintf(' Flux EDI (10^%g cm^{-1}s^{-1}sr^{-1})',log10(1/flux_scale));
    irf_legend(hca,sprintf('mms%g',ic),[0.02 0.98],[0 0 0])
    hold(hca,'on')
    plot(hca,[0 max([hca.XLim(2) hca.YLim(2)])],[0 max([hca.XLim(2) hca.YLim(2)])])
    hold(hca,'off')
    axis(hca,'square')
  end
end