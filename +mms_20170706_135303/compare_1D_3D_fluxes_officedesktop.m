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

%v_edi_edges = [v_edi_minus,v_edi_plus]*1e-3;

pa_edi_edges = -45:11.25:45;
pa_edi_centers = [(-45+11.25):11.25:0 0:11.25:(45-11.25)];
az_edges = 0:10:360;

% FPI energy and corresponding velocity, closest to EDI
E_fpi = double(ePDist1.elim(500).depend{1}(1,:)); % eV
v_fpi = sqrt(2*units.e*E_fpi./units.me); % m/s
dE_fpi = 65.01; % eV

E_fpi_plus = double(E_fpi + ePDist1.elim(500).ancillary.delta_energy_plus(1,:));
E_fpi_minus = double(E_fpi - ePDist1.elim(500).ancillary.delta_energy_minus(1,:));
v_fpi_plus = sqrt(2*units.e*E_fpi_plus./units.me); % m/s
v_fpi_minus = sqrt(2*units.e*E_fpi_minus./units.me); % m/s
v_fpi_plusminus = v_fpi_plus-v_fpi_minus;
dv_fpi_minus = v_fpi_minus - v_fpi;
dv_fpi_plus = v_fpi_plus - v_fpi;
dv_fpi = dv_fpi_plus - dv_fpi_minus; % m/s

v_fpi_edges = sort(-[v_fpi_minus,v_fpi_plus]);
v_edi_edges = sort(-[v_edi_minus,v_edi_plus]);

%%
tint_eh = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z');
tint_zoom = tint_zoom + [-1 1];

% PDist.reduce only takes in the center energy and define vg_edges from
% that. So to work around that we need to give two energies equally spaced
% from the desired energy, and half the desired bin width. Thereafter take
% the mean of the two bins.
c_eval('eDred?_ = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg'',sort(-[v_edi_minus v_edi_plus]+0.25*v_edi_plusminus*[1 -1])*1e-3);',1:4)
c_eval('eDred? = irf.ts_scalar(eDred?_.time,mean(eDred?_.data,2)); eDred?.units = eDred?_.units;',1:4)
c_eval('eDred?_edges_edi = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg_edges'',v_fpi_edges*1e-3);',1:4)
c_eval('eDred?_edges_fpi = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg_edges'',v_edi_edges*1e-3);',1:4)
%c_eval('eFlred? = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg'',([v_edi_minus v_edi_plus]+0.25*v_edi_plusminus*[1 -1])*1e-3);',1:4)
 dv = v_edi_plusminus; % m/s

c_eval('eFlux?_red_edges_edi = eDred?_edges_edi*v_edi*v_edi_plusminus; eFlux_red?_edges_edi = eFlux_red?_edges_edi*1e-4; eFlux_red?_edges_edi.units = ''1/s/m^2'';',1:4)
c_eval('eFlux?_red_edges_fpi = eDred?_edges_fpi*v_fpi*v_fpi_plusminus; eFlux_red?_edges_fpi = eFlux_red?_edges_fpi*1e-4; eFlux_red?_edges_fpi.units = ''1/s/m^2'';',1:4)
%c_eval('eFlux?_red_edges_edi = eDred?_edges_edi.flux; eFlux_red?_edges_edi = eFlux_red?_edges_edi*1e-4; eFlux_red?_edges_edi.units = ''1/s/m^2'';',1:4)
%c_eval('eFlux?_red_edges_fpi = eDred?_edges_fpi.flux; eFlux_red?_edges_fpi = eFlux_red?_edges_fpi*1e-4; eFlux_red?_edges_fpi.units = ''1/s/m^2'';',1:4)
c_eval('eFlux?_red = irf.ts_scalar(eDred?.time,sum(eDred?.data,2)*v_edi*dv); eFlux_red?.units = ''1/s/m^2'';',1:4)
c_eval('eFlux?_red = eDred?*v_edi*dv; eFlux_red?.units = ''1/s/m^2'';',1:4)
c_eval('eFlux?_red = eFlux_red?*1e-4; eFlux_red?.units = ''1/s/cm^2'';',1:4)
%%
c_eval('ePitch?_!bins = ePDist?.tlim(tint_zoom).pitchangles(dmpaB?,180 + [-11.25*! 0]);',1:4,1:4)
c_eval('ePitch?_4 = ePDist?.tlim(tint_zoom).pitchangles(dmpaB?,180 -11.25*[4 3 2 1 0]);',1:4)
c_eval('eFlux?_4 = ePitch?_4.flux;',1:4)
c_eval('eFlux?_4_Eedi = eFlux?_4.elim(E_edi);',1:4)
c_eval('ePitch?_16 = ePDist?.tlim(tint_zoom).pitchangles(dmpaB?,0:11.25:180);',1:4,1:4)
c_eval('eFlux?_fov_!bins = ePDist?.tlim(tint_zoom).flux.pitchangles(dmpaB?,180 + [-11.25*! 0]);',1:4,1:4)
c_eval('eFlux?_4_Eedi_ = irf.ts_scalar(eFlux?_4.time,squeeze(eFlux?_4.elim(E_edi).data));',1:4)

%%
npanels = 6;
nrows = 3;
ncols = 1;
[h,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical'); % horizontal

if 0 % phase space density fpi 4 bins cumsum
  hca = irf_panel('phase space density fpi 4 bins cumsum'); 
  irf_plot(hca,{ePitch1_4bins.elim(500),ePitch1_3bins.elim(500),ePitch1_2bins.elim(500),ePitch1_1bins.elim(500)},'comp')
  hca.YLabel.String = {'phase space density',sprintf('(%s)',ePitch1_1bins.units)};
  irf_legend(hca,{sprintf('[%.2f,%.2f]',ePitch1_4bins.ancillary.pitchangle_edges);...
                  sprintf('[%.2f,%.2f]',ePitch1_3bins.ancillary.pitchangle_edges);...
                  sprintf('[%.2f,%.2f]',ePitch1_2bins.ancillary.pitchangle_edges);...
                  sprintf('[%.2f,%.2f]',ePitch1_1bins.ancillary.pitchangle_edges)},[0.02 0.98])
end
if 0 % field-of-view flux fpi 4 bins cumsum
  hca = irf_panel('field-of-view flux fpi 4 bins cumsum');
  irf_plot(hca,{eFlux1_fov_4bins.elim(500),eFlux1_fov_3bins.elim(500),eFlux1_fov_2bins.elim(500),eFlux1_fov_1bins.elim(500)},'comp')
  hca.YLabel.String = {'field-of-view flux',sprintf('(%s)',eFlux_fov1_1bins.units)};
  irf_legend(hca,{sprintf('[%.2f,%.2f]',eFlux1_fov_4bins.ancillary.pitchangle_edges);...
                  sprintf('[%.2f,%.2f]',eFlux1_fov_3bins.ancillary.pitchangle_edges);...
                  sprintf('[%.2f,%.2f]',eFlux1_fov_2bins.ancillary.pitchangle_edges);...
                  sprintf('[%.2f,%.2f]',eFlux1_fov_1bins.ancillary.pitchangle_edges)},[0.02 0.98])
end
if 1
  hca = irf_panel('phase space density fpi'); 
  irf_plot(hca,{ePitch1_4.elim(E_edi)})
  hca.YLabel.String = {'phase space density',sprintf('(%s)',ePitch1_1bins.units)};
  irf_legend(hca,{sprintf('[%.2f,%.2f]',ePitch1_4.ancillary.pitchangle_edges(1:2));...
                  sprintf('[%.2f,%.2f]',ePitch1_4.ancillary.pitchangle_edges(2:3));...
                  sprintf('[%.2f,%.2f]',ePitch1_4.ancillary.pitchangle_edges(3:4));...
                  sprintf('[%.2f,%.2f]',ePitch1_4.ancillary.pitchangle_edges(4:5))},[0.02 0.98])
end
if 1
  hca = irf_panel('field-of-view flux fpi 4bins');
  irf_plot(hca,eFlux1_4.elim(E_edi))
  irf_legend(hca,{sprintf('[%.2f,%.2f]',eFlux1_4.ancillary.pitchangle_edges(1:2));...
                  sprintf('[%.2f,%.2f]',eFlux1_4.ancillary.pitchangle_edges(2:3));...
                  sprintf('[%.2f,%.2f]',eFlux1_4.ancillary.pitchangle_edges(3:4));...
                  sprintf('[%.2f,%.2f]',eFlux1_4.ancillary.pitchangle_edges(4:5))},[0.02 0.98])
end
if 0
  hca = irf_panel('field-of-view flux edi apar');
  irf_plot(hca,flux_135_180_mms1)
end
if 1
  hca = irf_panel('field-of-view flux edi apar');
  irf_plot(hca,ediFlux1.palim([90 180]))
end
if 1
  hca = irf_panel('field-of-view flux edi par');
  irf_plot(hca,ediFlux1.palim([0 90]))
end
if 0
  hca = irf_panel('field-of-view flux edi scaled');
  irf_plot(hca,flux_135_180_scaled_mms1)
end
if 0
  hca = irf_panel('reduced dist');
  irf_plot(hca,{eDred1,eDred2,eDred3,eDred4},'comp')
  hca.YLabel.String = {'reduced','phase space density',sprintf('(%s)',eDred1_.units)};
end
if 0
hca = irf_panel('reduced dist edi edges');
irf_plot(hca,{eDred1_edges_edi,eDred2_edges_edi,eDred3_edges_edi,eDred4_edges_edi},'comp')
hca.YLabel.String = {'reduced','phase space density',sprintf('(%s)',eDred1_.units)};
irf_legend(hca,{sprintf('(EDI bin) v = [%.0f,%.0f] km/s, dv = %.0f km/s',eDred1_edges_edi.ancillary.v_edges*1e-3,diff(eDred1_edges_edi.ancillary.v_edges*1e-3))},[0.02 0.98])

end
if 0
hca = irf_panel('reduced dist fpi edges');
irf_plot(hca,{eDred1_edges_fpi,eDred2_edges_fpi,eDred3_edges_fpi,eDred4_edges_fpi},'comp')
hca.YLabel.String = {'reduced','phase space density',sprintf('(%s)',eDred1_.units)};
irf_legend(hca,{sprintf('(closest FPI bin) v = [%.0f,%.0f] km/s, dv = %.0f km/s',eDred1_edges_fpi.ancillary.v_edges*1e-3,diff(eDred1_edges_fpi.ancillary.v_edges*1e-3))},[0.02 0.98])

end
if 0
hca = irf_panel('reduced flux');
irf_plot(hca,{eFlux_red1,eFlux_red2,eFlux_red3,eFlux_red4},'comp')
hca.YLabel.String = {'flux of reduced','distribution',sprintf('(%s)',eFlux_red1.units)};

end
if 0
hca = irf_panel('reduced flux edi edges');
irf_plot(hca,{eFlux1_red_edges_edi,eFlux2_red_edges_edi,eFlux3_red_edges_edi,eFlux4_red_edges_edi},'comp')
hca.YLabel.String = {'flux of reduced','distribution',sprintf('(%s)',eFlux_red1.units)};
irf_legend(hca,{sprintf('(EDI bin) v = [%.0f,%.0f] km/s, dv = %.0f km/s',eFlux1_red_edges_edi.ancillary.v_edges*1e-3,diff(eFlux1_red_edges_edi.ancillary.v_edges*1e-3))},[0.02 0.98])

end
if 0
hca = irf_panel('reduced flux fpi edges');
irf_plot(hca,{eFlux1_red_edges_fpi,eFlux2_red_edges_fpi,eFlux3_red_edges_fpi,eFlux4_red_edges_fpi},'comp')
hca.YLabel.String = {'flux of reduced','distribution',sprintf('(%s)',eFlux_red1.units)};
irf_legend(hca,{sprintf('(closest FPI bin) v = [%.0f,%.0f] km/s, dv = %.0f km/s',eFlux1_red_edges_fpi.ancillary.v_edges*1e-3,diff(eFlux1_red_edges_fpi.ancillary.v_edges*1e-3))},[0.02 0.98])
end           

irf_zoom(h,'x',tint_zoom)
irf_zoom(h,'y')

isub = 1;
if 1 % mean pa distribution during the interval, showing the pa bins and reduced velocity edges
  hca = h2(isub); isub = isub + 1;
  tint_pa = tint_eh+0.0*[-1 1];
  hpitch = ePitch1_16.plot_pad_polar(hca,'tint',tint_pa,'scpot',scPot1.resample(ePitch1_16));  
  hca.Title.String = sprintf('%s - %s',tint_pa(1).utc,tint_pa(2).utc);
  c_eval('irf_pl_mark(h(?),tint_pa)',1:numel(h))
  axis(hca,'square')
  hold(hca,'on')
  plot(hca,hca.XLim,[log10([E_fpi_plus;E_fpi_minus]); -log10([E_fpi_plus;E_fpi_minus])]*[1 1],'k')
  hold(hca,'off')
end
if 1 % mean pa distribution during the interval, showing the pa bins and reduced velocity edges
  hca = h2(isub); isub = isub + 1;
  tint_pa = tint_eh+0.0*[-1 1];
  hpitch = ePitch1_16.flux.plot_pad_polar(hca,'tint',tint_pa,'scpot',scPot1.resample(ePitch1_16));  
  hca.Title.String = sprintf('%s - %s',tint_pa(1).utc,tint_pa(2).utc);
  c_eval('irf_pl_mark(h(?),tint_pa)',1:numel(h))
  axis(hca,'square')
  hold(hca,'on')
  plot(hca,hca.XLim,[log10([E_fpi_plus;E_fpi_minus]); -log10([E_fpi_plus;E_fpi_minus])]*[1 1],'k')
  hold(hca,'off')
end
if 1 % mean pa distribution during the interval, showing the pa bins and reduced velocity edges
  hca = h2(isub); isub = isub + 1;
  tint_pa = tint_eh+0.0*[-1 1];
  hpitch = ePitch1_4.plot_pad_polar(hca,'tint',tint_pa,'scpot',scPot1.resample(ePitch1_16));  
  hca.Title.String = sprintf('%s - %s',tint_pa(1).utc,tint_pa(2).utc);
  c_eval('irf_pl_mark(h(?),tint_pa)',1:numel(h))
  axis(hca,'square')
  hold(hca,'on')
  plot(hca,hca.XLim,[log10([E_fpi_plus;E_fpi_minus]); -log10([E_fpi_plus;E_fpi_minus])]*[1 1],'k')
  hold(hca,'off')
end
if 0
  hca = h2(isub); isub = isub + 1;
  scatter(hca,eFlux1_fov_1bins.elim(500).data,eFlux_red1.data)
  if 1
    hold(hca,'on')
    scatter(hca,eFlux1_fov_2bins.elim(500).data,eFlux_red1.data)
    scatter(hca,eFlux1_fov_3bins.elim(500).data,eFlux_red1.data)
    scatter(hca,eFlux1_fov_4bins.elim(500).data,eFlux_red1.data)
    hold(hca,'off')
    irf_legend(hca,{sprintf('[%.2f,%.2f]',eFlux1_fov_1bins.ancillary.pitchangle_edges);...
                sprintf('[%.2f,%.2f]',eFlux1_fov_2bins.ancillary.pitchangle_edges);...
                sprintf('[%.2f,%.2f]',eFlux1_fov_3bins.ancillary.pitchangle_edges);...
                sprintf('[%.2f,%.2f]',eFlux1_fov_4bins.ancillary.pitchangle_edges)},[1.01 0.98])
  end

  hca.XLabel.String = {'field-of-view flux',sprintf('(%s)',eFlux_fov1_1bins.units)};
  hca.YLabel.String = {'flux of reduced distribution',sprintf('(%s)',eFlux_red1.units)};
  axis(hca,'equal','square')
  hca.XLim = [0 7]*1e6;
  %hca.YLim = hca.XLim;
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1
  
end