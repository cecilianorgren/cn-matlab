%load /Users/cecilia/Data/20170706_135303_basic_eh
units = irf_units;

E_fpi = 484;
v_fpi_484 = sqrt(2*units.e*E_fpi/units.me); % m/s


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


  
%% Figure
nrows = 2;
ncols = 3;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

f_flux_scale = 1e-4; % from s-1 m-2 > s-1 cm-2

for mms_id = 1:2 %mms_id = 1;
  
  
  c_eval('epdist = ePDist?.tlim(tint_phi);',mms_id)
  c_eval('flux0 = flux0_mms?.tlim(tint_phi);',mms_id)
  c_eval('flux180 = flux180_mms?.tlim(tint_phi);',mms_id)
  c_eval('ts_f_fpi = ef1D?.tlim(tint_phi);',mms_id)
  v_fpi = mean(ts_f_fpi.depend{1},1);
  f_fpi = ts_f_fpi.data;

  E_fpi_plus = epdist.ancillary.energy + epdist.ancillary.delta_energy_plus;
  E_fpi_minus = epdist.ancillary.energy - epdist.ancillary.delta_energy_minus;
  v_fpi_plus = sqrt(2*units.e*E_fpi_plus./units.me); % m/s
  v_fpi_minus = sqrt(2*units.e*E_fpi_minus./units.me); % m/s
  dv_fpi = v_fpi_plus - v_fpi_minus; % m/s
  dv_fpi_64 = [flipdim(dv_fpi,1) dv_fpi]; % m/s
  c_eval('v_fpi_ = ef1D?.tlim(tint_phi).depend{1};',mms_id) % km/s
  v_fpi_ = v_fpi_*1e3; % m/s 
  c_eval('ts_flux_fpi = ef1D?.tlim(tint_phi);',mms_id) 
  ts_flux_fpi.data = ts_flux_fpi.data.*abs(v_fpi_).*dv_fpi_64*f_flux_scale;
  
  c_eval('ts_flux_fpi_rescale = ef1D?.tlim(tint_phi);',mms_id)  
  ts_flux_fpi_rescale.data = ts_flux_fpi_rescale.data.*abs(v_fpi_).*dv_edi*f_flux_scale;
  
  flux_fpi = ts_flux_fpi.data;
  v_flux_fpi = mean(ts_f_fpi.depend{1},1);

  iE480 = find(abs(abs(v_fpi_(1,:))-v_fpi_480)==min(abs(abs(v_fpi_(1,:))-v_fpi_480)));
  ts_fpi_flux_edi = irf.ts_scalar(ts_flux_fpi.time,ts_flux_fpi.data(:,iE480(end:-1:1)));
  ts_fpi_flux_edi_rescale = irf.ts_scalar(ts_flux_fpi_rescale.time,ts_flux_fpi_rescale.data(:,iE480(end:-1:1)));
  
  %ts_fpi_flux_edi_rescale = ts_flux_fpi;
  %ts_fpi_flux_edi_rescale.data = ts_fpi_flux_edi_rescale.data.*abs(v_fpi_).*dv_edi*f_flux_scale;
  %ts_fpi_flux_edi = irf.ts_scalar(ts_flux_fpi.time,ts_flux_fpi.data(:,iE480(end:-1:1)));

  v_fpi = mean(ts_f_fpi.depend{1},1);
  f_fpi = ts_f_fpi.data;

  if 1 % FPI, 1D reduced 
    hca = h(isub); isub = isub + 1;
    v_scale = 1e-3;
    f_scale = 1e-8;
    hlines = plot(hca,v_fpi*v_scale,f_fpi*f_scale);
    hca.YLabel.String = {'f (s^1/m^4)'};
    hca.YLabel.String = {'f (s^1/cm^4)'};
    hca.XLabel.String = {'v (10^3 km/s)'};
    hca.XLim = [-40 40];
    fpi_utc = ts_f_fpi.time.utc;
    str_lines = {...
      sprintf('-- fpi: %s',fpi_utc(1,12:23));...
      sprintf('-- fpi: %s',fpi_utc(2,12:23));...
      sprintf('-- fpi: %s',fpi_utc(3,12:23));...
      sprintf('-- fpi: %s',fpi_utc(4,12:23))};  
    irf_legend(hca,str_lines,[0.99 0.99])  
    set(hca,'ColorOrder',zeros(10,3))  
    if 1 % EDI velocities                
      hold(hca,'on')
      all_edi_plusminus = [v_edi_minus;  v_edi_plus;...
                 -v_edi_minus; -v_edi_plus]*[1 1];     
       plot(hca,all_edi_plusminus*1e-6,hca.YLim,'k-.')
       irf_legend(hca,'EDI',[0.55 + 0.5*v_edi_plus*1e-6/hca.XLim(2) 0.5],[0 0 0])     
      hold(hca,'off')
    end
    hca.Title.String = sprintf('MMS %g',mms_id);
  end
  if 1 % FPI, 1D reduced, flux
    hca = h(isub); isub = isub + 1;
    v_scale = 1e-3;
    hlines = plot(hca,v_flux_fpi*v_scale,flux_fpi);
    hca.YLabel.String = {'f (s^1/m^4)'};
    hca.YLabel.String = {'flux (s^{-1}cm^{-2})'};
    hca.XLabel.String = {'v (10^3 km/s)'};
    hca.XLim = [-40 40];
    fpi_utc = ts_f_fpi.time.utc;
    str_lines = {...
      sprintf('-- fpi: %s',fpi_utc(1,12:23));...
      sprintf('-- fpi: %s',fpi_utc(2,12:23));...
      sprintf('-- fpi: %s',fpi_utc(3,12:23));...
      sprintf('-- fpi: %s',fpi_utc(4,12:23))};  
    irf_legend(hca,str_lines,[0.99 0.99])  
    set(hca,'ColorOrder',zeros(10,3))  
    if 1 % EDI velocities                
      hold(hca,'on')
      all_edi_plusminus = [v_edi_minus;  v_edi_plus;...
                 -v_edi_minus; -v_edi_plus]*[1 1];     
       plot(hca,all_edi_plusminus*1e-6,hca.YLim,'k-.')
       irf_legend(hca,'EDI',[0.55 + 0.5*v_edi_plus*1e-6/hca.XLim(2) 0.5],[0 0 0])     
      hold(hca,'off')
    end
    hca.Title.String = sprintf('MMS %g',mms_id);
  end
  if 1 % EDI flux at node 1
    hca = h(isub); isub = isub + 1; 
    %set(hca,'ColorOrder',mms_colors('matlab'))  
    irf_plot(hca,{flux0,flux180},'comp')
    hold(hca,'on')  
    %irf_plot(hca,ts_fpi_flux_edi,'o')
    irf_plot(hca,ts_fpi_flux_edi_rescale,'o')    
    hold(hca,'off')
    irf_legend(hca,{'edi 0^o','edi 180^o','fpi 0^o','fpi 180^o'},[0.98 0.99])
    irf_zoom(hca,'x',tint_phi)
    hca.YLabel.String = sprintf('flux (%s)',flux0.units);
    hca.YLabel.String = sprintf('flux (%s)','cm^{-2}s^{-1}');
    hca.YLabel.Interpreter = 'tex';
    hca.Title.String = sprintf('MMS %g',mms_id);    
  end
end

%%
it = 0;
if 0 % Flat skymap at energu Efpi (= 480)
  hca = h(isub); isub = isub + 1;
  it = it + 1;
  [ax,hcb] = mms.plot_skymap(hca,epdist(it),'energy',E_fpi,'flat','tint',epdist1.time(it));
end
if 0 % Flat skymap at energu Efpi (= 480)
  hca = h(isub); isub = isub + 1;
  it = it + 1;
  [ax,hcb] = mms.plot_skymap(hca,epdist1(it),'energy',E_fpi,'flat','tint',epdist1.time(it));
end
if 0 % Flat skymap at energu Efpi (= 480)
  hca = h(isub); isub = isub + 1;
  it = it + 1;
  [ax,hcb] = mms.plot_skymap(hca,epdist1(it),'energy',E_fpi,'flat','tint',epdist1.time(it));
end
if 0 % Flat skymap at energu Efpi (= 480)
  hca = h(isub); isub = isub + 1;
  it = it + 1;
  [ax,hcb] = mms.plot_skymap(hca,epdist1(it),'energy',E_fpi,'flat','tint',epdist1.time(it));
end


%% Plot with the steradian thing checked out
mms_id = 1;
nrows = 4;
ncols = 5;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

c_eval('ePDist = ePDist?.tlim(tint_phi);',mms_id)
c_eval('dmpaB = dmpaB?.resample(ePDist);',mms_id)
c_eval('scPot = scPot?.resample(ePDist);',mms_id)

eFlux = ePDist.vd3v('scpot',scPot,'sr');
eFluxPitch = eFlux.pitchangles(dmpaB,pitchangles_edges,'meanorsum','sum_weighted');

eFluxPitch15 = ePDist.vd3v('scpot',scPot,'sr').pitchangles(dmpaB,[0 15],'meanorsum','sum_weighted');


pitchangles_edges = [0 (11.25/2):((180-2*5.75)/15):(180-11.25/2) 180];
pitchangles_edges8 = [0 11.25:22.5:(180-11.25) 180];
pitchangles_edges8 = [0 15:30:165 180];

ePitch = ePDist.pitchangles(dmpaB,pitchangles_edges,'meanorsum','mean');
ePitch8 = ePDist.pitchangles(dmpaB,pitchangles_edges8,'meanorsum','mean');

eSr = ePitch.solidangle;
eSr8 = ePitch8.solidangle;

ePitchSum = ePDist.pitchangles(dmpaB,pitchangles_edges,'meanorsum','sum');
ePitchSum8 = ePDist.pitchangles(dmpaB,pitchangles_edges8,'meanorsum','sum');

ePitchSumWeight = ePDist.pitchangles(dmpaB,pitchangles_edges,'meanorsum','sum_weighted');
ePitchSumWeight8 = ePDist.pitchangles(dmpaB,pitchangles_edges8,'meanorsum','sum_weighted');

ePitchSumSr = ePitchSum; ePitchSumSr.data = ePitchSum.data./eSr.data; ePitchSumSr.units = sprintf('%s/sr',ePitchSumSr.units);
ePitchSumSr8 = ePitchSum8; ePitchSumSr8.data = ePitchSum8.data./eSr8.data; ePitchSumSr8.units = sprintf('%s/sr',ePitchSumSr8.units); 

ePitchSumWeightSr = ePitchSumWeight; ePitchSumWeightSr.data = ePitchSumWeight.data./eSr.data; ePitchSumWeightSr.units = sprintf('%s/sr',ePitchSumWeightSr.units);

ePitchMean = ePDist.pitchangles(dmpaB,pitchangles_edges,'meanorsum','mean');
ePitchMean8 = ePDist.pitchangles(dmpaB,pitchangles_edges8,'meanorsum','mean');

ePitchFlux = ePitch.pa_flux;
ePitchFluxSum = ePitchSum.pa_flux;
ePitchFluxSumWeight = ePitchSumWeight.pa_flux('sr');
ePitchFluxMean = ePitchMean.pa_flux;

ePitchFluxSr = ePitch.pa_flux('sr');
ePitchFluxSrSum = ePitchSum.pa_flux('sr');
ePitchFluxSrMean = ePitchMean.pa_flux('sr');

ePitchFlux8 = ePitch8.pa_flux;
ePitchFluxSum8 = ePitchSum8.pa_flux;
ePitchFluxSumWeight8 = ePitchSumWeight8.pa_flux('sr');
ePitchFluxMean8 = ePitchMean8.pa_flux;

ePitchFluxSr8 = ePitch8.pa_flux('sr');
ePitchFluxSrSum8 = ePitchSum8.pa_flux('sr');
ePitchFluxSrMean8 = ePitchMean8.pa_flux('sr');


eint = [000 40000];
vint = [-Inf Inf];
vint_pa = v_edi*sind(5.6250)*[-1 1]*1e-3;
c_eval('scPot = scPot?.resample(ePDist?);',mms_id)
lowerelim = irf.ts_scalar(ePDist.time,repmat(50,1,ePDist.length));
ePDist1D    = ePDist.elim(eint).reduce('1D',dmpaB.resample(ePDist).norm,'vint',vint   ,'scpot',scPot.resample(ePDist),'lowerelim',lowerelim);
ePDist1D_pa = ePDist.elim(eint).reduce('1D',dmpaB.resample(ePDist).norm,'vint',vint_pa,'scpot',scPot.resample(ePDist),'lowerelim',lowerelim);

it  = 1; % four times during tint_phi
if 0 % Flat skymap at energu E_fpi
  hca = h(isub); isub = isub + 1;
  [ax,hcb] = mms.plot_skymap(hca,ePDist(it),'energy',E_fpi,'flat','tint',ePDist.time(it),'vectors',{dmpaB.norm.data,'B'});    
end
if 0 % ePitch scpot corrected
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitch(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';  
  irf_legend(hca,{ePitch.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]); 
end
if 0 % ePitch scpot NOT corrected
  hca = h(isub); isub = isub + 1;
  ax = ePitch(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitch.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFlux scpot corrected
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchFlux(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';   
  irf_legend(hca,{ePitchFlux.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFlux scpot NOT corrected
  hca = h(isub); isub = isub + 1;
  ax = ePitchFlux(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFlux.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end

if 0 % ePitchMean scpot corrected pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchMean(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';  
  irf_legend(hca,{ePitchMean.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]); 
end
if 1 % ePitchMean scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchMean(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchMean.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFluxMean scpot corrected pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchFluxMean(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';   
  irf_legend(hca,{ePitchFluxMean.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchFluxMean scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchFluxMean(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxMean.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFluxMean scpot corrected pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchFluxSrMean(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';   
  irf_legend(hca,{ePitchFluxSrMean.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchFluxMean scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchFluxSrMean(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxSrMean.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % eSr scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = eSr(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';  
end

hca = h(isub); isub = isub + 1;

if 0 % ePitchSum scpot corrected pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchSum(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';  
  irf_legend(hca,{ePitchSum.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]); 
end
if 1 % ePitchSum scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchSum(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchSum.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchSumSr scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;  
  ePitchSumSr(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{sprintf('%s/sr',ePitchSumSr.ancillary.meanorsum)},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFluxSum scpot corrected pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchFluxSum(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';   
  irf_legend(hca,{ePitchFluxSum.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchFluxSum scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchFluxSum(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxSum.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFluxSum scpot corrected pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchFluxSrSum(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';   
  irf_legend(hca,{ePitchFluxSrSum.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchFluxSum scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchFluxSrSum(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxSrSum.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchFluxSum/eSr scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ePitchPlot = ePitchFluxSum;
  ePitchPlot.data = ePitchPlot.data./eSr.data;
  ePitchPlot(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{sprintf('%s/sr',ePitchPlot.ancillary.meanorsum)},[0.01 0.99],'color',[0 0 0]);
end

if 0 % ePitch scpot corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitch8(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';  
  irf_legend(hca,{ePitch8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]); 
end
if 1 % ePitch scpot NOT corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitch8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitch8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFlux scpot corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchFlux8(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';   
  irf_legend(hca,{ePitchFlux8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchFlux scpot NOT corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchFlux8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFlux8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end

if 0 % pitchangle flux per steradian
  hca = h(isub); isub = isub + 1;
  ax = loglog(hca,ePitch8(1).depend{1},squeeze(ePitch8(1).pa_flux('sr').data(1,:,[1 end]))');    
end
if 0 % ePitch scpot corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchMean8(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';  
  irf_legend(hca,{ePitchMean8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]); 
end
if 0 % ePitch scpot NOT corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchMean8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchMean8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFlux scpot corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchFluxMean8(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';   
  irf_legend(hca,{ePitchFluxMean8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFlux scpot NOT corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchFluxMean8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxMean8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitch scpot corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchSum8(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';  
  irf_legend(hca,{ePitchSum8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]); 
end
if 0 % ePitch scpot NOT corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchSum8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchSum8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFlux scpot corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  [ax,plot0_data] = ePitchFluxSum8(it).plot_pad_polar(hca,'scpot',scPot(it)); 
  hca.Title.String = 'Corrected for scPot';   
  irf_legend(hca,{ePitchFluxSum8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchFlux scpot NOT corrected 8 pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchFluxSum8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxSum8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end

if 0 % pitchangle psd
  hca = h(isub); isub = isub + 1;
  ax = loglog(hca,ePitch(1).depend{1},squeeze(ePitch(1).data(1,:,[1 end]))');
end
if 0 % pitchangle flux
  hca = h(isub); isub = isub + 1;
  ax = loglog(hca,ePitch(1).depend{1},squeeze(ePitch(1).pa_flux.data(1,:,[1 end]))');    
end
if 0 % pitchangle flux per steradian
  hca = h(isub); isub = isub + 1;
  ax = loglog(hca,ePitch(1).depend{1},squeeze(ePitch(1).pa_flux('sr').data(1,:,[1 end]))');    
end
if 0 % pitchangle flux
  hca = h(isub); isub = isub + 1;
  ax = loglog(hca,ePitch8(1).depend{1},squeeze(ePitch8(1).pa_flux.data(1,:,[1 end]))');    
end

%% Comparing mean and sum/sr
nrows = 3;
ncols = 4;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;


if 1 % ePitchMean scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchMean(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchMean.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchSum scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchSum(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchSum.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchSumWeight scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchSumWeight(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchSumWeight.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % eSr scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = eSr(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';  
  irf_legend(hca,{'sr'},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchSumSr scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;  
  ePitchSumSr(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{sprintf('%s/sr',ePitchSumSr.ancillary.meanorsum)},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchSumSr scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;  
  ePitchSumWeightSr(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{sprintf('%s/sr',ePitchSumWeightSr.ancillary.meanorsum)},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchFluxSum scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchFluxSumWeight(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxSumWeight.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end

if 1 % ePitchMean8 scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchMean8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchMean8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchSum8 scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchSum8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchSum8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchSumWeight scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchSumWeight8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchSumWeight8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % eSr8 scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = eSr8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';  
  irf_legend(hca,{'sr'},[0.01 0.99],'color',[0 0 0]);
end
if 1 % ePitchFluxSum scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ax = ePitchFluxSumWeight8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxSumWeight8.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end
if 0 % ePitchSumSr8 scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;  
  ePitchSumSr8(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{sprintf('%s/sr',ePitchSumSr8.ancillary.meanorsum)},[0.01 0.99],'color',[0 0 0]);
end

if 1 % compare 180 pa to reduced flux
  hca = h(isub); isub = isub + 1;
  ax = plot(hca,ePitchFluxSumWeight(it).depend{1}, squeeze(ePitchFluxSumWeight(it).data(1,:,[end]))',...
                ePitchFluxSumWeight8(it).depend{1},squeeze(ePitchFluxSumWeight8(it).data(1,:,[end]))'...
    );    
  hca.XLabel.String = 'E (eV)';
end
if 1 % compare 180 pa to reduced flux
  hca = h(isub); isub = isub + 1;
  ax = plot(hca,ePDist1D(it).depend{1}, squeeze(ePDist1D(it).data(1,:))',...
                  ePDist1D_pa(it).depend{1},squeeze(ePDist1D_pa(it).data(1,:))'...
    );    
  hca.XLabel.String = 'v (km/s)';
end


if 1 % eFluxPitch
  hca = h(isub); isub = isub + 1;
  ax = eFluxPitch(it).plot_pad_polar(hca);
  hca.Title.String = 'Flux corrected for scPot, not plot';
  irf_legend(hca,{eFluxPitch.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0]);
end

%% Comparing pitch angle fluxes to reduced fluxes
mms_id = 1;
nrows = 3;
ncols = 3;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

c_eval('ePDist = ePDist?.tlim(tint_phi);',mms_id)
%ePDist_nobg = mms.remove_edist_background(ePDist);
ePDist = ePDist_nobg;
c_eval('dmpaB = dmpaB?.resample(ePDist);',mms_id)
c_eval('scPot = scPot?.resample(ePDist);',mms_id)
c_eval('flux180 = flux180_mms?.tlim(tint_phi).resample(ePDist);',mms_id)

pitchangles_edges = [0 (11.25/2):((180-2*5.75)/15):(180-11.25/2) 180];
pitchangles_edges8 = [0 11.25:22.5:(180-11.25) 180];
pitchangles_edges8 = [0 15:30:165 180];

thetas = 5:5:35;


if 1 % dont use scpot to correct for stuff  
  eFlux = ePDist.vd3v('sr');
  eFluxPitch = eFlux.pitchangles(dmpaB,pitchangles_edges,'meanorsum','sum_weighted');
  c_eval('eFluxPitch_0_? = ePDist.vd3v(''sr'').pitchangles(dmpaB,[180-? 180],''meanorsum'',''sum_weighted'');',thetas)  

  ePitch = ePDist.pitchangles(dmpaB,pitchangles_edges,'meanorsum','sum_weighted');
  ePitchFlux = ePitch.pa_flux();
  ePitchFluxSr = ePitch.pa_flux('sr');

  eint = [000 40000];
  vint = [-Inf Inf];
  c_eval('scPot = scPot?.resample(ePDist?);',mms_id)
  lowerelim = irf.ts_scalar(ePDist.time,repmat(00,1,ePDist.length));
  ePDist1D = ePDist.elim(eint).reduce('1D',dmpaB.resample(ePDist).norm,'vint',vint);
  c_eval('ePDist1D_0_? = ePDist.elim(eint).reduce(''1D'',dmpaB.resample(ePDist).norm,''vint'',v_edi*sind(?)*[-1 1]*1e-3);',thetas)
  
  eFlux1D = ePDist1D.red_flux;
  c_eval('eFlux1D_0_? = ePDist1D_0_?.red_flux;',thetas);
else
  eFlux = ePDist.vd3v('scpot',scPot,'sr');
  eFluxPitch = eFlux.pitchangles(dmpaB,pitchangles_edges,'meanorsum','sum_weighted');
  c_eval('eFluxPitch_0_? = ePDist.vd3v(''scpot'',scPot,''sr'').pitchangles(dmpaB,[180-? 180],''meanorsum'',''sum_weighted'');',thetas)

  ePitch = ePDist.pitchangles(dmpaB,pitchangles_edges,'meanorsum','sum_weighted');
  ePitchFlux = ePitch.pa_flux('scpot',scPot);
  ePitchFluxSr = ePitch.pa_flux('scpot',scPot,'sr');

  eint = [000 40000];
  vint = [-Inf Inf];  
  c_eval('scPot = scPot?.resample(ePDist?);',mms_id)
  lowerelim = irf.ts_scalar(ePDist.time,repmat(00,1,ePDist.length));
  ePDist1D      = ePDist.elim(eint).reduce('1D',dmpaB.resample(ePDist).norm,'vint',vint                      ,'scpot',scPot.resample(ePDist),'lowerelim',lowerelim);
  c_eval('ePDist1D_0_? = ePDist.elim(eint).reduce(''1D'',dmpaB.resample(ePDist).norm,''vint'',v_edi*sind(?)*[-1 1]*1e-3,''scpot'',scPot.resample(ePDist),''lowerelim'',lowerelim);',thetas)
  
  eFlux1D = ePDist1D.red_flux;
  c_eval('eFlux1D_0_? = ePDist1D_0_?.red_flux;',thetas);
end

% This is not how it should be
% 
eFlux1D_allthetas_data = nan(eFlux1D.length,size(eFlux1D.depend{1},2),numel(thetas));
for ith = 1:numel(thetas)
  c_eval('eFlux1D_allthetas_data(:,:,ith) = eFlux1D_0_?.data;',thetas(ith))
end
v_edi_ind = find(abs(eFlux1D_allthetas.userData.depend{1}(1,:)--v_edi*1e-3)==min(abs(eFlux1D_allthetas.userData.depend{1}(1,:)--v_edi*1e-3)));
eFlux1D_allthetas = irf.ts_scalar(eFlux1D.time,eFlux1D_allthetas_data);
eFlux1D_allthetas.userData.depend{1} = eFlux1D.depend{1};
eFlux1D_allthetas.userData.depend{2} = thetas;
eFlux1D_allthetas.units = eFlux1D.units;

eFlux1D_allthetas_vedi = irf.ts_scalar(eFlux1D_allthetas.time,squeeze(eFlux1D_allthetas.data(:,v_edi_ind,:)));
eFlux1D_allthetas_vedi.userData.depend = {eFlux1D_allthetas.userData.depend{2}};
eFlux1D_allthetas_vedi.units = eFlux1D.units;
  

it  = 1; % four times during tint_phi
if 1 % ePitch scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ePitch(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitch.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0],'interpreter','none');
end
if 1 % ePitchFlux scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ePitchFluxSr(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxSr.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0],'interpreter','none');
end
if 1 % eFluxPitch scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  eFluxPitch(it).plot_pad_polar(hca);
  hca.Title.String = 'Flux corrected for scPot, plot no';
  irf_legend(hca,{eFluxPitch.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0],'interpreter','none');
end
if 0 % ePitchFluxSr-eFluxPitch
  hca = h(isub); isub = isub + 1;
  flux_diff = ePitchFluxSr-eFluxPitch;
  flux_diff.data = abs(flux_diff.data);
  flux_diff(it).plot_pad_polar(hca);
  hca.Title.String = 'Flux corrected for scPot, plot no';
  irf_legend(hca,{flux_diff.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0],'interpreter','none');
end
if 1 % ePitchFlux scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ePitchFluxSr(it).plot_pad_polar(hca);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxSr.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0],'interpreter','none');
  if 1 % vint used for fred
    hold(hca,'on')
    colors = mms_colors('matlab');
    ii = 0;
    for theta_int = 5:5:35
      ii = ii + 1;
      plot_int = log10(E_edi)*sind(theta_int);
      lines_pos = plot(hca,(plot_int*[1 1]),hca.YLim,'linewidth',1.5);
      lines_neg = plot(hca,-(plot_int*[1 1]),hca.YLim,'linewidth',1.5);
      lines_pos.Color = colors(ii,:);
      lines_neg.Color = colors(ii,:);
    end      
    hold(hca,'off')
  end
  if 1 % solid angle used for pitch angle fluxes
    hold(hca,'on')
    colors = mms_colors('matlab');    
    for ith = 1:numel(thetas)
      theta_int = thetas(ith);      
      plotx = log10(E_edi)*sind(theta_int);
      ploty = log10(E_edi)*cosd(theta_int);
      
      lines_sr1 = plot(hca,[0  plotx],[0 -ploty],'linewidth',1.5);
      lines_sr2 = plot(hca,[0 -plotx],[0 -ploty],'linewidth',1.5);      
      lines_sr1.Color = colors(ith,:);
      lines_sr2.Color = colors(ith,:);            
    end      
    hold(hca,'off')
  end
end
if 1 % ePitchFlux scpot NOT corrected pa's
  hca = h(isub); isub = isub + 1;
  ePitchFluxSr(it).plot_pad_polar(hca,'scpot',scPot);
  hca.Title.String = 'Not corrected for scPot';
  irf_legend(hca,{ePitchFluxSr.ancillary.meanorsum},[0.01 0.99],'color',[0 0 0],'interpreter','none');
  if 1 % vint used for fred
    hold(hca,'on')
    colors = mms_colors('matlab');
    ii = 0;
    for theta_int = 5:5:35
      ii = ii + 1;
      plot_int = log10(E_edi-scPot(it).data)*sind(theta_int);
      lines_pos = plot(hca,(plot_int*[1 1]),hca.YLim,'linewidth',1.5);
      lines_neg = plot(hca,-(plot_int*[1 1]),hca.YLim,'linewidth',1.5);
      lines_pos.Color = colors(ii,:);
      lines_neg.Color = colors(ii,:);
    end      
    hold(hca,'off')
  end
  if 1 % solid angle used for pitch angle fluxes
    hold(hca,'on')
    colors = mms_colors('matlab');    
    for ith = 1:numel(thetas)
      theta_int = thetas(ith);      
      plotx = log10(E_edi-scPot(it).data)*sind(theta_int);
      ploty = log10(E_edi-scPot(it).data)*cosd(theta_int);
      
      lines_sr1 = plot(hca,[0  plotx],[0 -ploty],'linewidth',1.5);
      lines_sr2 = plot(hca,[0 -plotx],[0 -ploty],'linewidth',1.5);      
      lines_sr1.Color = colors(ith,:);
      lines_sr2.Color = colors(ith,:);            
    end      
    hold(hca,'off')
  end
end 
if 1 % Reduced f: for a few different vints
  hca = h(isub); isub = isub + 1;
  plot(hca,ePDist1D_0_5(it).depend{1}, squeeze(ePDist1D_0_5(it).data)',...
           ePDist1D_0_10(it).depend{1}, squeeze(ePDist1D_0_10(it).data)',...
           ePDist1D_0_15(it).depend{1}, squeeze(ePDist1D_0_15(it).data)',...
           ePDist1D_0_20(it).depend{1}, squeeze(ePDist1D_0_20(it).data)',...
           ePDist1D_0_25(it).depend{1}, squeeze(ePDist1D_0_25(it).data)',...
           ePDist1D_0_30(it).depend{1}, squeeze(ePDist1D_0_30(it).data)',...
           ePDist1D_0_35(it).depend{1}, squeeze(ePDist1D_0_35(it).data)',...
           'Marker','*'...
      );    
  legend(hca,{'\pm 5^o','\pm 10^o','\pm 15^o','\pm 20^o','\pm 25^o','\pm 30^o','\pm 35^o'},'location','best')
  hca.XLabel.String = sprintf('v (%s)',ePDist1D_0_05.ancillary.vint_units);
  hca.YLabel.String = sprintf('f(v) (%s)',ePDist1D_0_05.units);
  hca.Title.String = 'Reduced distribution';
  hca.XLim = 0.3e5*[-1 1];
end
if 1 % Reduced: Flux/sr for a few different vints
  hca = h(isub); isub = isub + 1;
  plot(hca,eFlux1D_0_5(it).depend{1}, abs(squeeze(eFlux1D_0_5(it).data))',...
           eFlux1D_0_10(it).depend{1}, abs(squeeze(eFlux1D_0_10(it).data))',...
           eFlux1D_0_15(it).depend{1}, abs(squeeze(eFlux1D_0_15(it).data))',...
           eFlux1D_0_20(it).depend{1}, abs(squeeze(eFlux1D_0_20(it).data))',...
           eFlux1D_0_25(it).depend{1}, abs(squeeze(eFlux1D_0_25(it).data))',...
           eFlux1D_0_30(it).depend{1}, abs(squeeze(eFlux1D_0_30(it).data))',...
           eFlux1D_0_35(it).depend{1}, abs(squeeze(eFlux1D_0_35(it).data))',...
           'Marker','*'...
      );    
  if 1 % add downsampled edi flux
    hold(hca,'on')
    plot(hca,hca.XLim,flux180.data(it)*1e4*[1 1],'k--')
    hold(hca,'off')
  end  
  if 1 % add edi energy
    hold(hca,'on')
    loglog(hca,-v_edi*1e-3*[1 1],hca.YLim,'k--',v_edi*1e-3*[1 1],hca.YLim,'k--')
    hold(hca,'off')
  end
  legend(hca,{'\theta = \pm 5^o','\pm 10^o','\pm 15^o','\pm 20^o','\pm 25^o','\pm 30^o','\pm 35^o'},'location','best')
  hca.XLabel.String = 'v (.)';
  hca.YLabel.String = sprintf('Flux (%s)',eFlux1D_0_05.units);
  hca.Title.String = {'Flux (parallel) of reduced distribution','for different vints = v_{edi}sind(\theta)'};
  hca.XLim = 0.3e5*[-1 1];
  hca.YLim = 3e10*[0 1];
end
if 1 % Log Pitch angle: Flux for a few different solid angles
  hca = h(isub); isub = isub + 1;
  loglog(hca,eFluxPitch_0_5(it).depend{1}, squeeze(eFluxPitch_0_5(it).data)',...
             eFluxPitch_0_10(it).depend{1}, squeeze(eFluxPitch_0_10(it).data)',...
             eFluxPitch_0_15(it).depend{1}, squeeze(eFluxPitch_0_15(it).data)',...
             eFluxPitch_0_20(it).depend{1}, squeeze(eFluxPitch_0_20(it).data)',...
             eFluxPitch_0_25(it).depend{1}, squeeze(eFluxPitch_0_25(it).data)',...
             eFluxPitch_0_30(it).depend{1}, squeeze(eFluxPitch_0_30(it).data)',...
             eFluxPitch_0_35(it).depend{1}, squeeze(eFluxPitch_0_35(it).data)',...
             'Marker','*'...
    );   
  if 1 % add downsampled edi flux
    hold(hca,'on')
    loglog(hca,hca.XLim,flux180.data(it)*[1 1],'k--')
    hold(hca,'off')
  end
  if 1 % add edi energy
    hold(hca,'on')
    loglog(hca,E_edi*[1 1],hca.YLim,'k--')
    hold(hca,'off')
  end
  
  legend(hca,{'\pm 5^o','\pm 10^o','\pm 15^o','\pm 20^o','\pm 25^o','\pm 30^o','\pm 35^o'},'location','best')
  hca.XLabel.String = 'E (eV)';
  hca.YLabel.String = sprintf('Flux (%s)',eFluxPitch_0_05.units);
  hca.Title.String = 'Antiparallel flux for different solid angles';
end
if 1 % Lin Pitch angle: Flux for a few different solid angles
  hca = h(isub); isub = isub + 1;
  plot(hca,eFluxPitch_0_5(it).depend{1}, squeeze(eFluxPitch_0_5(it).data)',...
             eFluxPitch_0_10(it).depend{1}, squeeze(eFluxPitch_0_10(it).data)',...
             eFluxPitch_0_15(it).depend{1}, squeeze(eFluxPitch_0_15(it).data)',...
             eFluxPitch_0_20(it).depend{1}, squeeze(eFluxPitch_0_20(it).data)',...
             eFluxPitch_0_25(it).depend{1}, squeeze(eFluxPitch_0_25(it).data)',...
             eFluxPitch_0_30(it).depend{1}, squeeze(eFluxPitch_0_30(it).data)',...
             eFluxPitch_0_35(it).depend{1}, squeeze(eFluxPitch_0_35(it).data)',...
             'Marker','*'...
    );   
  if 1 % add downsampled edi flux
    hold(hca,'on')
    plot(hca,hca.XLim,flux180.data(it)*[1 1],'k--')
    hold(hca,'off')
  end
  if 1 % add edi energy
    hold(hca,'on')
    plot(hca,E_edi*[1 1],hca.YLim,'k--')
    hold(hca,'off')
  end
  hca.XScale = 'log';
  legend(hca,{'\pm 5^o','\pm 10^o','\pm 15^o','\pm 20^o','\pm 25^o','\pm 30^o','\pm 35^o'},'location','best')
  hca.XLabel.String = 'E (eV)';
  hca.YLabel.String = sprintf('Flux (%s)',eFluxPitch_0_05.units);
  hca.Title.String = 'Antiparallel flux for different solid angles';
end

