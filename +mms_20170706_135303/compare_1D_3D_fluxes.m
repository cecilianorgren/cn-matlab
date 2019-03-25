units = irf_units;
tint_zoom = irf.tint('2017-07-06T13:54:00.000Z/2017-07-06T13:54:10.000Z');
tint_zoom = irf.tint('2017-07-06T13:53:55.000Z/2017-07-06T13:54:15.000Z');
% EDI energy and corresponding velocity
E_edi = 500; % eV
v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
dE_edi = 25; % eV

E_edi_plus = E_edi + dE_edi;
E_edi_minus = E_edi - dE_edi;
v_edi_plus = sqrt(2*units.e*E_edi_plus./units.me); % m/s
v_edi_minus = sqrt(2*units.e*E_edi_minus./units.me); % m/s
v_edi_plusminus = v_edi_plus-v_edi_minus;
dv_edi_minus = v_edi_minus - v_edi;
dv_edi_plus = v_edi_plus - v_edi;
dv_edi = dv_edi_plus - dv_edi_minus; % m/s

v_edi_edges = [v_edi_minus,v_edi_plus]*1e-3;
pa_edi_edges = -45:11.25:45;
pa_edi_centers = [(-45+11.25):11.25:0 0:11.25:(45-11.25)];
az_edges = 0:10:360;

tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z');
tint_zoom = tint_zoom + [-1 1];

c_eval('ePitch?_flux_edi = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',1:4)
c_eval('ePitch?_fpi = ePDist?.pitchangles(dmpaB?,16);',1:4);
c_eval('ePitch?_flux_fpi = ePitch?_fpi.flux;',1:4);
c_eval('ePitch?_flux_fpi_bin_closest_to_edi = ePitch?_flux_fpi.elim(500);',1:4);
%c_eval('ePitch?_flux_fpi_edi_range = ePitch?_flux_fpi.rebin({});',1:4);

%c_eval('eRed?_fpi_edi_range = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg_edges'',(-[v_edi_minus v_edi_plus])*1e-3);',1:4)
c_eval('eRed?_fpi_edi_range_180 = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg_edges'',([-v_edi_plus -v_edi_minus ])*1e-3);',1)
c_eval('eRed?_fpi_flux_edi_range = eRed?_fpi_edi_range.flux_red;',1)

%% Plot
npanels = 4;
h = irf_plot(npanels);
isub = 1;

legends_edi = {};
for ipa=1:8 legends_edi{ipa} = sprintf('[%.2f %.2f]',ePitch1_flux_edi.depend{2}(ipa)-ePitch1_flux_edi.ancillary.delta_pitchangle_minus(ipa)*0.5,ePitch1_flux_edi.depend{2}(ipa)+ePitch1_flux_edi.ancillary.delta_pitchangle_plus(ipa)*0.5); end

if 0 % edi flux par
  hca = h(isub); isub = isub + 1;
  palim = [0 45];
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6})
  irf_legend(hca,legends_edi(1:4),[0.01 0.98])
  hca.YLabel.String = {'Flux edi','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.Interpreter = 'tex';
end
if 0 % edi fpi 0
  hca = h(isub); isub = isub + 1;
  palim = 0;
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch1_flux_fpi.palim(palim).elim(500)*1e-6},'comp')
  irf_legend(hca,{'EDI','FPI'},[0.01 0.98])
  hca.YLabel.String = {'Flux','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.Interpreter = 'tex';
end
if 0 % edi flux apar
  hca = h(isub); isub = isub + 1;
  palim = [135 180];
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6})
  irf_legend(hca,legends_edi(5:8),[0.01 0.98])
  hca.YLabel.String = {'Flux edi','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.Interpreter = 'tex';
end
if 1 % edi fpi 180
  hca = h(isub); isub = isub + 1;
  palim = 180;
  set(hca,'ColorOrder',pic_colors('xy'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch1_flux_fpi.palim(palim).elim(500)*1e-6},'comp')
  irf_legend(hca,{'EDI','FPI'},[0.01 0.98])
  hca.YLabel.String = {'Flux','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.Interpreter = 'tex';
end
if 1 % edi fpi 180, edi resampled to fpi
  hca = h(isub); isub = isub + 1;
  palim = 180;
  set(hca,'ColorOrder',pic_colors('xy'))
  irf_plot(hca,{ePitch1_flux_edi.resample(ePitch1_flux_fpi).palim(palim)*1e-6,ePitch1_flux_fpi.palim(palim).elim(500)*1e-6},'comp')
  irf_legend(hca,{'EDI','FPI'},[0.01 0.98])
  hca.YLabel.String = {'Flux','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.Interpreter = 'tex';
end
if 1 % fpi red edi e range
  hca = h(isub); isub = isub + 1;
  palim = 180;
  set(hca,'ColorOrder',pic_colors('z'))
  irf_plot(hca,{eRed1_fpi_edi_range_180.flux_red*1e-6},'comp')
  irf_legend(hca,{sprintf('v = [%.0f, %.0f] km/s',eRed1_fpi_edi_range_180.ancillary.v_edges(1)*1e-3,eRed1_fpi_edi_range_180.ancillary.v_edges(2)*1e-3)},[0.01 0.05])
  hca.YLabel.String = {'Reduced flux','FPI','(10^6 s^{-1}cm^{-2})'};
  hca.YLabel.Interpreter = 'tex';
end
if 1 % fpi red edi e range, edi fpi 180
  hca = h(isub); isub = isub + 1;
  palim = 180;
  set(hca,'ColorOrder',pic_colors('xyz'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch1_flux_fpi.palim(palim).elim(500)*1e-6,eRed1_fpi_edi_range_180.flux_red*1e-6},'comp')
  irf_legend(hca,{'EDI fov',...
    'FPI fov',...
    sprintf('FPI red: v = [%.0f, %.0f] km/s',eRed1_fpi_edi_range_180.ancillary.v_edges(1)*1e-3,eRed1_fpi_edi_range_180.ancillary.v_edges(2)*1e-3)},...
    [0.01 0.98])
  hca.YLabel.String = {'Flux','red (10^6 s^{-1}cm^{-2})','fov (10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.Interpreter = 'tex';
end
if 0 % fpi red edi e range apar, edi fpi 0, just compare for debugging
  hca = h(isub); isub = isub + 1;
  palim = 0;
  set(hca,'ColorOrder',pic_colors('xyz'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch1_flux_fpi.palim(palim).elim(500)*1e-6,eRed1_fpi_edi_range_180.flux_red*1e-6},'comp')
  irf_legend(hca,{'EDI fov',...
    'FPI fov',...
    sprintf('FPI red: v = [%.0f, %.0f] km/s',eRed1_fpi_edi_range_180.ancillary.v_edges(1)*1e-3,eRed1_fpi_edi_range_180.ancillary.v_edges(2)*1e-3)},...
    [0.01 0.98])
  hca.YLabel.String = {'Flux','red (10^6 s^{-1}cm^{-2})','fov (10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.Interpreter = 'tex';
end
irf_zoom(h,'x',tint)

%%
%c_eval('eDred? = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg'',([v_edi_minus v_edi_plus]+0.25*v_edi_plusminus*[1 -1])*1e-3);',1:4)

%c_eval('eFlred? = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg'',([v_edi_minus v_edi_plus]+0.25*v_edi_plusminus*[1 -1])*1e-3);',1:4)
%dv = v_edi_plusminus; % m/s
%c_eval('eDred?_ = irf.ts_scalar(eDred?.time,sum(eDred?.data,2)); eDred?_.units = eDred?.units;',1:4)
%c_eval('eFlux_red? = irf.ts_scalar(eDred?.time,sum(eDred?.data,2)*v_edi*dv); eFlux_red?.units = ''1/s/m^2'';',1:4)
%c_eval('eFlux_red? = eFlux_red?*1e-4; eFlux_red?.units = ''1/s/cm^2'';',1:4)

%c_eval('ePitch?_!bins = ePDist?.tlim(tint_zoom).pitchangles(dmpaB?,180 + [-11.25*! 0]);',1:4,1:4)
%c_eval('eFlux_fov?_!bins = ePDist?.tlim(tint_zoom).flux.pitchangles(dmpaB?,180 + [-11.25*! 0]);',1:4,1:4)

%%
npanels = 4;
[h,h2] = initialize_combined_plot(npanels,1,1,0.6,'vertical'); % horizontal

hca = irf_panel('phase space density');
irf_plot(hca,{ePitch1_1bins.elim(500),ePitch1_2bins.elim(500),ePitch1_3bins.elim(500),ePitch1_4bins.elim(500)},'comp')
hca.YLabel.String = {'phase space density',sprintf('(%s)',ePitch1_1bins.units)};


hca = irf_panel('field-of-view flux');
irf_plot(hca,{eFlux_fov1_1bins.elim(500),eFlux_fov1_2bins.elim(500),eFlux_fov1_3bins.elim(500),eFlux_fov1_4bins.elim(500)},'comp')
hca.YLabel.String = {'field-of-view flux',sprintf('(%s)',eFlux_fov1_1bins.units)};

hca = irf_panel('reduced dist');
irf_plot(hca,{eDred1_,eDred2_,eDred3_,eDred4_},'comp')
hca.YLabel.String = {'reduced','phase space density',sprintf('(%s)',eDred1_.units)};

hca = irf_panel('reduced flux');
irf_plot(hca,{eFlux_red1,eFlux_red2,eFlux_red3,eFlux_red4},'comp')
hca.YLabel.String = {'flux of reduced distribution',sprintf('(%s)',eFlux_red1.units)};

irf_zoom(h,'x',tint_zoom)

isub = 1;
hca = h2(isub); isub = isub + 1;
scatter(hca,eFlux_fov1_1bins.elim(500).data,eFlux_red1.data)
if 0
  hold(hca,'on')
  scatter(hca,eFlux_fov1_2bins.elim(500).data,eFlux_red1.data)
  scatter(hca,eFlux_fov1_3bins.elim(500).data,eFlux_red1.data)
  scatter(hca,eFlux_fov1_4bins.elim(500).data,eFlux_red1.data)
  hold(hca,'off')
end

hca.XLabel.String = {'field-of-view flux',sprintf('(%s)',eFlux_fov1_1bins.units)};
hca.YLabel.String = {'flux of reduced distribution',sprintf('(%s)',eFlux_red1.units)};
axis(hca,'equal','square')
hca.XLim = [0 3]*1e6;
hca.YLim = [0 3]*1e6;
hca.Box = 'on';
hca.XGrid = 'on';
hca.YGrid = 'on';