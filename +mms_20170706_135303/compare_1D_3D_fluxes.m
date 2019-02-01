units = irf_units;

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

c_eval('eDred? = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg'',([v_edi_minus v_edi_plus]+0.25*v_edi_plusminus*[1 -1])*1e-3);',1:4)
%c_eval('eFlred? = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg'',([v_edi_minus v_edi_plus]+0.25*v_edi_plusminus*[1 -1])*1e-3);',1:4)
dv = v_edi_plusminus; % m/s
c_eval('eDred?_ = irf.ts_scalar(eDred?.time,sum(eDred?.data,2)); eDred?_.units = eDred?.units;',1:4)
c_eval('eFlux_red? = irf.ts_scalar(eDred?.time,sum(eDred?.data,2)*v_edi*dv); eFlux_red?.units = ''1/s/m^2'';',1:4)
c_eval('eFlux_red? = eFlux_red?*1e-4; eFlux_red?.units = ''1/s/cm^2'';',1:4)

c_eval('ePitch?_!bins = ePDist?.tlim(tint_zoom).pitchangles(dmpaB?,180 + [-11.25*! 0]);',1:4,1:4)
c_eval('eFlux_fov?_!bins = ePDist?.tlim(tint_zoom).flux.pitchangles(dmpaB?,180 + [-11.25*! 0]);',1:4,1:4)

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