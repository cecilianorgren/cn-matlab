mms_id = 1;
units = irf_units;
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint_fred = irf.tint('2017-07-06T13:54:05.55Z/2017-07-06T13:54:05.65Z');

%% EDI ranges
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

%% Load data
% Magnetic field
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',mms_id);

% Spacecraft potential
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',mms_id);
%c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);

% Skymap distributions
c_eval('tic; [iPDist?,iPDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',mms_id)
c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',mms_id)

%% Make reduced distribution
eint = [000 Inf];
vint = [-Inf Inf];
vg_cart = (-100:2:100)*1e3;
vg_par = vg_cart;
lowerelim = 40;
nMC = 200;

c_eval('edist = ePDist?.tlim(tint_fred);',mms_id)
c_eval('scpot = scPot?.resample(edist);',mms_id)

c_eval('ePara = dmpaB?.resample(edist).norm;',mms_id)
ePerp1 = ePara.cross(irf.ts_vec_xyz(ePara.time,repmat([1 0 0],ePara.length,1))).norm;
ePerp2 = ePara.cross(ePerp1).norm;

tic; ef1D = edist.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg_par,'nMC',nMC); toc % reduced distribution along B
tic; ef2D_parperp1 = edist.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg_cart,'base','cart','nMC',nMC); toc 
tic; ef2D_parperp2 = edist.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg_cart,'base','cart','nMC',nMC); toc
tic; ef2D_perp1perp2 = edist.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg_cart,'base','cart','nMC',nMC); toc

tic; ef1D_edi = edist.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg_edges',[v_edi_minus v_edi_plus]*1e-3,'nMC',nMC); toc % reduced distribution along B
tic; ef1D_edi_vint = edist.reduce('1D',ePara,'vint',v_edi_plus*[-1 1]*1e-3*sind(2*11.25),'scpot',scpot,'lowerelim',lowerelim,'vg_edges',[v_edi_minus v_edi_plus]*1e-3,'nMC',nMC); toc % reduced distribution along B

ePitch = edist.pitchangles(ePara,16);
ePitchPar = edist.pitchangles(ePara,[168.75 180]);

%% Plot
for it = 2 % plot for some given time index
  time_str = edist(it).time.utc;
  figure(21)
  nrows = 3;
  ncols = 3;
  npanels = ncols*nrows;
  for ipanel = 1:npanels
    h(ipanel) = subplot(nrows,ncols,ipanel);
  end
  isub = 1;
  
  if 1 % par - timseseries 
    hca = h(isub); isub = isub + 1;
    irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
    irf_timeaxis(hca)
    %hca.XLabel.String = 'v_{||} (10^3 km/s)';
    %hca.YLabel.String = 'f (s/m^4)';
  end
  if 1 % par
    hca = h(isub); isub = isub + 1;
    plot(hca,ef1D(it).depend{1}*1e-3,ef1D(it).data);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'f (s/m^4)';
    hca.Title.String = time_str(1:23);
  end
  if 1 % parperp1
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ef2D_parperp1(it).plot_plane(hca);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  end
  if 1 % par1perp2
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ef2D_parperp2(it).plot_plane(hca,'contourf',[]);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
  end
  if 1 % perp1perp2
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ef2D_perp1perp2(it).plot_plane(hca,'contour',5);
    hca.XLabel.String = 'v_{\perp1} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp2} (10^3 km/s)';
  end
  if 1 % pitch angle spectrogram
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ePitch(it).plot_pad_polar(hca,'scpot',scpot);
    hca.XLabel.String = 'v_{\perp} (10^3 km/s)';
    hca.YLabel.String = 'v_{||} (10^3 km/s)';
    axis(hca,'square');
  end
  if 1 % pitch angle spectrogram - par
    hca = h(isub); isub = isub + 1;
    [h_surf,h_axis,h_all] = ePitchPar(it).plot_pad_polar(hca,'scpot',scpot);
    hca.XLabel.String = 'v_{\perp} (10^3 km/s)';
    hca.YLabel.String = 'v_{||} (10^3 km/s)';
    axis(hca,'square');
  end
    
  for ipanel = [2]
    h(ipanel).CLim = [-13 -9.5];
    axis(h(ipanel),'square');
    h(ipanel).XLim = vg_par([1 end])*1e-3*0.5;
    %h(ipanel).YLim = [0 4e-3];
  end  
  for ipanel = [3:5]
    h(ipanel).CLim = [-13 -9.5];
    axis(h(ipanel),'square');
    h(ipanel).XLim = vg_cart([1 end])*1e-3*1;
    h(ipanel).YLim = vg_cart([1 end])*1e-3*1;
  end  
end