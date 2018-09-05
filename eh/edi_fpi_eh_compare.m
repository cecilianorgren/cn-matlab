%load /Users/cecilia/Data/20170706_135303_basic_eh
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

% Figur
nrows = 2;
ncols = 3;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

Efpi = 480;
v_fpi_480 = sqrt(2*units.e*Efpi/units.me); % m/s
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
  [ax,hcb] = mms.plot_skymap(hca,epdist(it),'energy',Efpi,'flat','tint',epdist1.time(it));
end
if 0 % Flat skymap at energu Efpi (= 480)
  hca = h(isub); isub = isub + 1;
  it = it + 1;
  [ax,hcb] = mms.plot_skymap(hca,epdist1(it),'energy',Efpi,'flat','tint',epdist1.time(it));
end
if 0 % Flat skymap at energu Efpi (= 480)
  hca = h(isub); isub = isub + 1;
  it = it + 1;
  [ax,hcb] = mms.plot_skymap(hca,epdist1(it),'energy',Efpi,'flat','tint',epdist1.time(it));
end
if 0 % Flat skymap at energu Efpi (= 480)
  hca = h(isub); isub = isub + 1;
  it = it + 1;
  [ax,hcb] = mms.plot_skymap(hca,epdist1(it),'energy',Efpi,'flat','tint',epdist1.time(it));
end
