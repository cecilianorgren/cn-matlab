%% Load datastore, save directory
localuser = datastore('local','user');
saveDirRoot = ['/Users/' localuser '/GoogleDrive/Research/cold_ion_pickup/'];
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   

%% Set spacecraft and tint
ic = 1;
tint = irf.tint('2017-07-06T15:46:23.00Z/2017-07-06T15:48:53.00Z');
tint = irf.tint('2017-07-03T21:54:03.00Z/2017-07-03T21:55:53.00Z');
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');

% Event string
fileName = ePDist1.userData.GlobalAttributes.Logical_file_id;
fileNameSplit = strsplit(fileName{1},'_'); numName = fileNameSplit{6};
strEvent = sprintf('%s-%s-%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));

%% Load data
units = irf_units;
% Magnetic field
disp('Loading magnetic field...')
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);

% Electric field
disp('Loading electric field...')
c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
c_eval('E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint);',ic);

% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint);',ic);

% Skymap distributions
disp('Loading skymaps...')
%c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
%c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)

c_eval('tic; [iPDist?,iPDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)
c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)

% Remove all one-count "noise"
c_eval('iPDist?.data(iPDist?.data<iPDistErr?.data*1.1) = 0;',ic)
%c_eval('ePDist?.data(ePDist?.data<ePDistErr?.data*1.1) = 0;',ic)

% Pressure and temperature
disp('Loading pressure and temperature...');
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic);

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
disp('Loading density...');
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);

% Velocity
disp('Loading bulk velocities...');
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic)

disp('Done loading data.');

%% Prepare data
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s
c_eval('EiExB? = units.mp*gseVExB?.abs.resample(gseVi?.time).^2*1e6/2/units.eV;',ic)
c_eval('EeExB? = units.me*gseVExB?.abs.resample(gseVe?.time).^2*1e6/2/units.eV;',ic)
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('[gseVi?par,gseVi?perp] = irf_dec_parperp(gseB?,gseVi?); gseVi?par.name = ''Vi par''; gseVi?perp.name = ''Vi perp'';',ic)
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

% Make reduced distribution
eint = [000 40000];
vint = [-Inf Inf];

iDist = iPDist1.elim(eint);

%tic; ef1D = eDist.reduce('1D',eLine,'vint',vint,'scpot',scpot_lim); toc % reduced distribution along B
tic; if1Dpar = iDist.reduce('1D',dmpaB1.resample(iDist).norm); toc % reduced distribution along B
tic; if1Dx = iDist.reduce('1D',[1 0 0]); toc % reduced distribution along B
tic; if1Dy = iDist.reduce('1D',[0 1 0]); toc % reduced distribution along B
% there is an error if direction is exactly [0 0 1]. Due to cross with [0 0 1] to make a system orthogonal 
tic; if1Dz = iDist.reduce('1D',[0 1 100]); toc % reduced distribution along B

%lineVe = ve.dot(eLine); % projection of Vi on B
%lineVi = vi.dot(iLine); % projection of Vi on B

%% Plot
tintZoom = tint;
ic = 1;
npanels = 9;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
cmap = [1 1 1; pic_colors('waterfall')];
clim_ifred = [-3 -1];
doExB = 1;
doV = 1;

if 0 % B abs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B abs');
  irf_plot(hca,gseB1.abs);
  hca.YLabel.String = 'B (nT)';
end
if 0 % B GSM
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
  hca.YLabel.String = {'B^{GSM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 1 % B GSE
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 1 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  %irf_legend(hca,{'n_e','n_i'},[0.98 9.98])
  hca.YLabel.String = {'n_e','(cm^{-3})'};
end
if 0 % iDEF omni
  isub = isub + 1;
  hca = irf_panel('iDEF');  
  [hout,hcb] = irf_spectrogram(hca,iDist.deflux.omni.specrec,'log');
  if 0 % vExB energy
    hold(hca,'on')
    lineScpot = irf_plot(hca,EiExB1,'k');
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 0.5;
    hold(hca,'off')
  end
  colormap(hca,cmap)
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YLim(1) = 100;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 1 % iDEF omni, with ExB energy
  isub = isub + 1;
  hca = irf_panel('iDEF+EExB');  
  [hout,hcb] = irf_spectrogram(hca,iDist.deflux.omni.specrec,'log');
  if 1 % vExB energy
    hold(hca,'on')
    lineExB = irf_plot(hca,EiExB1,'k');
    lineExB.Marker = '.';
    lineExB.LineStyle = 'none';
    lineExB.MarkerSize = 1;
    lineExB.Color = [0 0 0]; %lineExB.LineWidth = 0.5;
    hold(hca,'off')
  end
  if 1 % spacecraft potential
    hold(hca,'on')
    lineScpot = irf_plot(hca,scPot1,'k');
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 0.5;
    hold(hca,'off')
  end
  colormap(hca,cmap)
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YLim(1) = 100;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLabel.Interpreter = 'tex';
end
if 0 % iPDist pa 64
  isub = isub + 1;
  hca = irf_panel('i PA e64 deflux lowe');  
  eint = [100 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
end
if 1 % i psd vpar
  isub = isub + 1;
  hca = irf_panel('iLine');
  irf_plot(hca,if1Dpar.specrec('velocity_1D'),'nocolorbar');
  hold(hca,'on')
  irf_plot(hca,{gseVi1par},'comp')
  %irf_plot(hca,gseVi1)
  colormap(hca,cmap)
  hold(hca,'off')
  hca.YLim = if1Dpar.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_{i,||}','(km/s)'}; 
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.CLim = clim_ifred;
end
if 1 % i psd vx
  isub = isub + 1;
  hca = irf_panel('iLine x');
  irf_plot(hca,if1Dx.specrec('velocity_1D'),'nocolorbar');
  if doV
    hold(hca,'on')
    irf_plot(hca,{gseVi1.x},'comp')
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    irf_plot(hca,{gseVExB1.resample(gseVi1).x},'comp')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLim = if1Dx.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_{i,x}','(km/s)'}; 
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.CLim = clim_ifred;
end
if 1 % i psd vy
  isub = isub + 1;
  hca = irf_panel('iLine y');
  irf_plot(hca,if1Dy.specrec('velocity_1D'),'nocolorbar');
  if doV
    hold(hca,'on')
    irf_plot(hca,{gseVi1.y},'comp')
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    irf_plot(hca,{gseVExB1.resample(gseVi1).y},'comp')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLim = if1Dy.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_{i,y}','(km/s)'}; 
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.CLim = clim_ifred;
end
if 1 % i psd vz
  isub = isub + 1;
  hca = irf_panel('iLine z');
  irf_plot(hca,if1Dz.specrec('velocity_1D'),'nocolorbar');
  if doV
    hold(hca,'on')
    irf_plot(hca,{gseVi1.z},'comp')
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    irf_plot(hca,{gseVExB1.resample(gseVi1).z},'comp')
    hold(hca,'off')
  end
    %irf_plot(hca,gseVi1)
  colormap(hca,cmap)
  hca.YLim = if1Dz.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_{i,z}','(km/s)'};  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.CLim = clim_ifred;
end
if 1 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % eDEF omni
  isub = isub + 1;
  hca = irf_panel('eDEF');
  [hout,hcb] = irf_spectrogram(hca,eDist.omni.specrec,'log');
  hold(hca,'on')
  lineScpot = irf_plot(hca,scpot_lim,'k');
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_e','(eV)'};   
  irf_legend(hca,[num2str(scpot_margin) 'V_{sc}'],[0.99 0.1],'color',0*[1 1 1])
end
if 0 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('eLine');
  irf_plot(hca,ef1D.specrec('velocity_1D'));
  hold(hca,'on')
  irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_e (km/s)'; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 0 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end

if 0 % Te par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
end
if 0 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % E perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x,gseE?perp.y,gseE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % E par, perp, downsampled to vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp par downsampled');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.resample(gseVi1).x,gseE?perp.resample(gseVi1).y,gseE?perp.resample(gseVi1).z,gseE?par.resample(gseVi1)},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\perp,x','\perp,y','\perp,z','||'},[0.98 0.9],'fontsize',12);  
end

irf_zoom(h,'x',tintZoom)
%irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h=irf_plot({gseB1,gseVi1,iPDist1.deflux.omni.specrec('energy'),f1D.specrec('velocity_1D')}); h(3).YScale = 'log'; %h(4).YLim = [-1000 1000];

