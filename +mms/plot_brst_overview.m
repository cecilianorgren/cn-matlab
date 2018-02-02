%% Load data
ic = 2:3;

tint = irf.tint('2015-11-12T06:42:04.00Z/2015-11-12T06:43:00.00Z'); %20151112064204
tint = irf.tint('2016-12-25T16:00:14.00Z',30); %20161225160014, or 20161225160144. Dec 25 2016 - at about 16.00
% time is redefined to brst whole interval
tint = irf.tint('2016-11-28T08:25:40.00Z/2016-11-28T08:27:30.00Z');
tint = irf.tint('2016-11-25T08:03:10.00Z/2016-11-25T08:03:50.00Z');
tint = irf.tint('2016-12-29T11:12:14.00Z/2016-12-29T11:13:53.00Z'); % mms1_sp_20161229_111214_111353
tint = irf.tint('2017-01-01T13:59:30.00Z/2017-01-01T14:00:00.00Z');
tint = irf.tint('2015-12-02T01:14:30.00Z/2015-12-02T01:15:15.00Z');
tint = irf.tint('2015-10-18T15:11:30.00Z/2015-10-18T15:12:45.00Z'); 

tint = irf.tint('2017-07-06T15:36:43.00Z/2017-07-06T15:41:13.00Z'); % 20170706153643
tint = irf.tint('2017-07-06T15:39:12.00Z/2017-07-06T15:40:13.00Z'); % 20170706153913
%tint = irf.tint('2017-07-06T15:41:33.00Z/2017-07-06T15:42:33.00Z'); % 20170706154133
%tint = irf.tint('2017-07-06T15:44:03.00Z/2017-07-06T15:45:33.00Z'); % 20170706154403
%tint = irf.tint('2017-07-06T15:46:23.00Z/2017-07-06T15:47:43.00Z'); % 20170706154623
%tint = irf.tint('2017-07-06T15:48:53.00Z/2017-07-06T15:50:43.00Z'); % 20170706154853
%tint = irf.tint('2017-07-06T15:51:13.00Z/2017-07-06T15:52:43.00Z'); % 20170706155113
%tint = irf.tint('2017-07-06T15:53:43.00Z/2017-07-06T15:54:43.00Z'); % 20170706155343
%tint = irf.tint('2017-07-06T15:56:03.00Z/2017-07-06T15:57:43.00Z'); % 20170706155603
%
%
%tint = irf.tint('2017-07-11T22:29:23.00Z/2017-07-11T22:31:43.00Z'); % 20170711222923
%tint = irf.tint('2017-07-11T22:33:23.00Z/2017-07-11T22:34:43.00Z'); % 20170711223323, major X-line moving tailward
%tint = irf.tint('2017-07-11T22:37:23.00Z/2017-07-11T22:38:43.00Z'); % 20170711223723

%tint = irf.tint('2017-07-09T10:45:13.00Z/2017-07-09T10:46:43.00Z'); % 20170709104513
%tint = irf.tint('2017-07-09T10:47:03.00Z/2017-07-09T10:48:43.00Z'); % 20170709104703
%tint = irf.tint('2017-07-09T10:49:03.00Z/2017-07-09T10:51:03.00Z'); % 20170709104903

% 
tint = irf.tint('2017-06-11T17:24:53.00Z',1); % 20170611172453
tint = irf.tint('2017-06-11T17:32:53.00Z',1); % 20170611173253
tint = irf.tint('2017-06-11T17:35:03.00Z',1); % 20170611173503
tint = irf.tint('2017-06-11T17:37:03.00Z',1); % 20170611173703
tint = irf.tint('2017-06-11T17:39:13.00Z',1); % 20170611173913
tint = irf.tint('2017-06-11T17:41:23.00Z',1); % 20170611174123
tint = irf.tint('2017-06-11T17:43:23.00Z',1); % 20170611174323
% tint = irf.tint('2017-06-11T17:52:23.00Z',1); % 20170611175223
% 
% tint = irf.tint('2017-06-19T10:41:23.00Z',1); % 20170619104123
% tint = irf.tint('2017-06-19T10:43:43.00Z',1); % 20170619104343
% tint = irf.tint('2017-06-19T10:46:13.00Z',1); % 20170619104613
% tint = irf.tint('2017-06-19T10:48:43.00Z',1); % 20170619104843
% tint = irf.tint('2017-06-19T10:51:13.00Z',1); % 20170619105113
% tint = irf.tint('2017-06-19T10:38:53.00Z',1); % 20170619103853 - not done
% tint = irf.tint('2017-06-19T10:36:23.00Z',1); % 20170619103623 - not done

% tint = irf.tint('2017-07-03T21:52:13.00Z',1); % 20170703215213
% tint = irf.tint('2017-07-03T21:54:03.00Z',1); % 20170703215403
% tint = irf.tint('2017-07-03T21:55:13.00Z',1); % 20170703215553
% tint = irf.tint('2017-07-03T21:57:43.00Z',1); % 20170703215743
% tint = irf.tint('2017-07-03T21:59:33.00Z',1); % 20170703215933
% tint = irf.tint('2017-07-03T22:01:23.00Z',1); % 20170703220123
% tint = irf.tint('2017-07-03T22:03:23.00Z',1); % 20170703220323
% tint = irf.tint('2017-07-03T22:05:33.00Z',1); % 20170703220533

tint = irf.tint('2015-10-16T13:06:00.00Z',60); % 20170703220533
tint = irf.tint('2017-06-17T20:24:00.00Z',60);

% Set datastore
mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/data/mms');
db_info = datastore('mms_db');   

% Particle distributions: electrons and ions
disp('Loading particle distributions...')
c_eval('tic; ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
c_eval('tic; iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)

% Setting tint from ePDist1
tint = ePDist1([1 ePDist1.length]).time;

% Make event directory
fileName = ePDist1.userData.GlobalAttributes.Logical_file_id;
fileNameSplit = strsplit(fileName{1},'_'); numName = fileNameSplit{6};
dirName = sprintf('%s-%s-%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
%eventPath = ['/Users/Cecilia/Research/Events/' dirName '/'];
eventPath = ['/Users/cno062/Research/Events/' dirName '/']; 
mkdir(eventPath)
    
% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB?scm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint); toc',ic);

% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);

% Load spacecraft position
disp('Loading spacecraft position...')
mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/data/mms');

try
R = mms.get_data('R_gse',tint);
if ~all([isfield(R,'gseR1') isfield(R,'gseR1') isfield(R,'gseR2') isfield(R,'gseR3') isfield(R,'gseR4')])  
  % not the right data fiels, try to load from irfu database instead
  db_info = datastore('mms_db');   
  mms.db_init('local_file_db','/data/mms');
  R = mms.get_data('R_gse',tint);
  mms.db_init('local_file_db',db_info.local_file_db_root);
end
if size(R.gseR1,2) == 4
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
else
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end
end
% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);

% Pressure and temperature
disp('Loading pressure and temperature...')
c_eval('tic; gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?); toc;',ic) 
c_eval('tic; gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?); toc;',ic) 
c_eval('tic; gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?); toc;',ic)

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
disp('Loading density...')
c_eval('tic; ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?); toc;',ic);
c_eval('tic; ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?); toc;',ic);

% Velocity
disp('Loading bulk velocities...')
c_eval('tic; gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)

% Current: use electron density to calculate current
units = irf_units;
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);

% Pitchangle distribution
c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,13); ePitch? = ePitch?.convertto(''s^3/km^6'');',ic)

% Feeps data
%c_eval('dobj = dataobj(''/Volumes/Nexus/data/mms?/feeps/brst/l2/electron/2017/07/11/mms?_feeps_brst_l2_electron_20170711222923_v5.5.1.cdf'')',ic)
%c_eval('feeps? = mms.variable2ts(get_variable(dobj,''mms?_epd_feeps_brst_l2_electron_spin''));',ic)
%mms1_epd_feeps_brst_l2_electron_spin

c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.units = ''mV/m'';',ic)
c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.units = ''mV/m'';',ic)

c_eval('gseEVexB? = gseE?.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''E+VexB'';',ic)
%% Plot overview figure with focus on electrons
npanels = 10;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
cmap = colormap('jet');

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 0 % ne ni
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 1 % J  
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint),gseVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x_{\perp}','y_{\perp}','z_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 1 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 0 % e DEF omni 32
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 1 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [100 30000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
try
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
%irf_zoom(h([1:3 6 10]),'y')
%hca = irf_panel('Te'); irf_zoom(hca,'y')
%hca = irf_panel('e PA e64 deflux lowe');  hca.YLim = [0 180];
%hca = irf_panel('n'); hca.YLim = [0 11];
%hca = irf_panel('J zoom'); hca.YLim = [-500 800];
%hca = irf_panel('Ve'); hca.YLim = [-900 500];
%hca = irf_panel('Te'); hca.YLim = [20 110]; hca.YTick = [20:20:120];
%hca = irf_panel('gradPe'); hca.YLim = [-2.2 2.2];
%hca = irf_panel('B brst'); hca.YLabel.String = {'B','(nT)'}
%hca = irf_panel('Ti'); hca.YLim = [300 700];
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
if 0
hmark = irf_pl_mark(h(1:6),tintBCS, 'yellow');
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end
for ii = 1:npanels;
  h(ii).FontSize = 12;
end
%irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))

%% Plot overview figure with focus on electrons, including single time electron distributions
npanels = 10;
tintZoom = tint;%irf.tint('2016-12-25T16:00:53.00Z/2016-12-25T16:00:58.00Z');
cmap = 'jet';
[h1,h2] = initialize_combined_plot(npanels,2,2,0.4,'vertical');
ic = 1;
iisub = 0;
cmap = colormap('jet');

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 1 % J  
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint),gseVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x_{\perp}','y_{\perp}','z_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 1 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 1 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [100 30000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h1(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h1,'x',tintZoom)
%irf_zoom(h([1:3 6 10]),'y')
%hca = irf_panel('Te'); irf_zoom(hca,'y')
%hca = irf_panel('e PA e64 deflux lowe');  hca.YLim = [0 180];
%hca = irf_panel('n'); hca.YLim = [0 11];
%hca = irf_panel('J zoom'); hca.YLim = [-500 800];
%hca = irf_panel('Ve'); hca.YLim = [-900 500];
%hca = irf_panel('Te'); hca.YLim = [20 110]; hca.YTick = [20:20:120];
%hca = irf_panel('gradPe'); hca.YLim = [-2.2 2.2];
%hca = irf_panel('B brst'); hca.YLabel.String = {'B','(nT)'}
%hca = irf_panel('Ti'); hca.YLim = [300 700];
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);
if 0
hmark = irf_pl_mark(h1(1:6),tintBCS, 'yellow');
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end
for ii = 1:npanels;
  h1(ii).FontSize = 12;
end
%h1(7).CLim = [7 8.5];
%% Plot single time particle distributions, 1 sc, 4 projections,
  
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('ePitch = ePitch?;',ic)

% Plot format input
vlim = 50*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [-3 1.5];  
  

c_eval('times = ePDist?.time;',ic)
tind = 786;%:1:920;953; 911;
tind = 1160;
tind = 1240;
tind = 1451:1500;%1350:1450;
for it = 1251:2:1381;%2260:2300%:1710%:2370;1420:1530;2340:2370;tind;589;%:650;550:650;%1240;tind;
  time = times(it);
  it
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,time.epochUnix','green');
  
  c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
  c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
  c_eval('hatExB = cross(hatE,hatB);',ic)
  par = hatB;
  perp1 = hatExB;
  perp2 = cross(par,perp1);  
  

  timeUTC = time.utc;      
  isub = 1;

  %vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B'};
  vectors = {[1 0 0],'x'; [0 1 0],'y'; [0 0 1],'z'};
  
  % Perpendicular plane    
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{perp2}','v_{||}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
  hca.Title.String = timeUTC(1:23);

  % B plane 1
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;par;-perp2]; vlabels = {'v_{ExB}','v_{||}','-v_{perp2}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % B plane 2
  hca = h2(isub); isub = isub + 1;
  xyz = [perp2;par;perp1]; vlabels = {'v_{perp2}','v_{||}','v_{ExB}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % Pitchangle distribution
  hca = h2(isub); isub = isub + 1;
  plot(hca,ePitch.depend{1}(it,:),squeeze(ePitch.data(it,:,[1 ceil(numel(ePitch.depend{2})/2) numel(ePitch.depend{2})])));
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' ePitch.units ')'];
  hca.XLabel.String = 'E (eV)';
  hca.XLim = [10 30000];
  legend(hca,{'0','90','180'})
  hca.YTick = 10.^[-4:5];
  hca.YLim = [1e-4 1e5];

  % ExB plane, with markings
  if 0
    vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z'};
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    hca.Title.String = '';
  end
  for ii = 1:4
    colormap(h2(ii),strCMap)
  end
  cn.print(['e_proj_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','path',eventPath)
end
