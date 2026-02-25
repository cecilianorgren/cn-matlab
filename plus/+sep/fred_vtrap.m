ic = 1;

% Set datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/data/mms');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
%% tint



tint = irf.tint('2017-07-18T01:40:23.00Z/2017-07-18T01:42:13.00Z');
tint_fred = irf.tint('2017-07-18T01:40:30.00Z',5);
tint_zoom = tint_fred;
%% Load data
% Particle distributions: electrons and ions
disp('Loading particle distributions...')
c_eval('tic; ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
c_eval('tic; iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)

% Setting tint from ePDist1
c_eval('tint = ePDist?([1 ePDist?.length]).time;',ic)

% Make event directory
c_eval('fileName = ePDist?.userData.GlobalAttributes.Logical_file_id;',ic)
fileNameSplit = strsplit(fileName{1},'_'); numName = fileNameSplit{6};
eventName = sprintf('%s-%s-%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
eventPath = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/']; 
mkdir(eventPath)
    
% Magnetic field
disp('Loading magnetic field...')
c_eval('dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',1:4);
c_eval('gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('gseB?scm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',1:4);

% Electric field
disp('Loading electric field...')
c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);

% Load spacecraft position
disp('Loading spacecraft position...')

try
  R = mms.get_data('R_gse',tint);
  if ~all([isfield(R,'gseR1') isfield(R,'gseR1') isfield(R,'gseR2') isfield(R,'gseR3') isfield(R,'gseR4')])  
    % not the right data fielsd, try to load from irfu database instead
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
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

% Pressure and temperature
disp('Loading pressure and temperature...')
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic)

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

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)

% Current, curlometer
if all([~isempty(gseB1) ~isempty(gseB2) ~isempty(gseB3) ~isempty(gseB4)])
  c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
  [Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
  gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
  gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
  gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';
else
  gseJcurl = irf.ts_vec_xyz(gseB1.time,gseB1.data*NaN);
end

%% Event path
fileName = ePDist1.userData.GlobalAttributes.Logical_file_id;
fileNameSplit = strsplit(fileName{1},'_'); numName = fileNameSplit{6};
dirName = sprintf('%s-%s-%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
dirNameMatlab = sprintf('+mms_%s%s%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
eventPath = ['/Users/' localuser '/Research/Events/' dirName '/'];
matlabPath = ['/Users/' localuser '/MATLAB/cn-matlab/' dirNameMatlab '/'];

%% Fred
eDist = ePDist1.tlim(tint_fred);
ve = gseVe1.tlim(eDist.time).resample(eDist);
 % keep in mind that this also affects the velocity at lower energies
scpot = scPot1.resample(eDist);
lowerelim = 50;
energies = sqrt(2*eDist.depend{1}(1,:)*units.eV/units.me)/1000; % km/s
vgmax = 70000;
%vg = [-energies(end:-1:1) 0 energies];
vg = -vgmax:1000:vgmax;
vg(abs(vg)>70000) = [];
tic; ef1D = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'scpot',scpot,'lowerelim',lowerelim,'vg',vg); toc % reduced distribution along B
disp('Done.')

%% Get wave properties
% Load data from txt-file
fid = fopen([matlabPath 'esw_properties.txt'],'r');
data_format_read = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f'; 
esw_data = textscan(fid,[data_format_read]);
fclose(fid)

all_t_center = EpochTT(char(esw_data{5}));
[all_t_center_sorted,ind_sorted] = all_t_center.sort;

for icell = 1:numel(esw_data)
  esw_data{icell} = esw_data{icell}(ind_sorted);
  esw_data{icell}(8) = [];
end

fid = fopen([matlabPath 'esw_properties_redo.txt'],'r');
  
tsVph = irf.ts_vec_xyz(char(esw_data{5}),[esw_data{6} esw_data{7} esw_data{8}]);
tsVphpar = tsVph.dot(gseB1.norm.resample(tsVph));
tsPhi = irf.ts_scalar(char(esw_data{5}),[esw_data{10} esw_data{11} esw_data{12} esw_data{13}]);
tsVtrap = irf.ts_scalar(char(esw_data{5}),sqrt(2*units.e*esw_data{10}/units.me)*1e-3);
c_eval('tsVtrap? = irf.ts_scalar(char(esw_data{5}),sqrt(2*units.e*esw_data{9+?}/units.me)*1e-3);',1:4)

%% Plot, including f proj and v phi and vtrap, for AGU talk 2018
ic = 1;
npanels = 4;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
tint_zoom = tint_zoom;

vmin = tsVphpar-tsVtrap;
vmax = tsVphpar+tsVtrap;
  
fmin = 10^(-7);

if 0 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('eLine');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %set(hca,'ColorOrder',mms_colors('122'))
  %irf_plot(hca,{tsVphpar,vmin,vmax},'comp');
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  %hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_e (km/s)'; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 0 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseVe?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};
  irf_legend(hca,sprintf('mms%g',ic),[0.98 0.9],'fontsize',12)
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'}; 
  irf_legend(hca,sprintf('mms%g',ic),[0.98 0.9],'fontsize',12)
end
if 1 % Phi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{tsPhi});
  c_eval('hh(?).Color = mms_colors(''?'')',1:4)
  hca.YLabel.String = {'\Phi_{||}','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[0.98 0.9],'fontsize',12);
end
hca = irf_panel('delete');
if 1 % v phase + trap
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phase velocity');
  
  ef_plot = ef1D;
  ef_plot.data(ef_plot.data<fmin) = NaN;
  [~,hcb] = irf_spectrogram(hca,ef_plot.specrec('velocity_1D','10^3 km/s'));
  hold(hca,'on')
  c_eval('vmin? = tsVphpar-tsVtrap?;',1:4)
  c_eval('vmax? = tsVphpar+tsVtrap?;',1:4)
  
  %set(hca,'ColorOrder',mms_colors('111223344'))
  set(hca,'ColorOrder',mms_colors('111111111'))
  vscale = 1e-3;
  %irf_plot(hca,{tsVphpar*vscale,vmin1*vscale,vmax1*vscale,vmin2*vscale,vmax2*vscale,vmin3*vscale,vmax3*vscale,vmin4*vscale,vmax4*vscale},'comp');
  hp = irf_plot(hca,{tsVphpar*vscale,vmin*vscale,vmax*vscale},'comp');
  
  %irf_patch(hca,{vmin,vmax})
  %hca.YLim = sort(real([max([vmax1.data; vmax2.data; vmax3.data; vmax4.data]) min([vmin1.data; vmin2.data; vmin3.data; vmin4.data])]));
  hold(hca,'off')
  hca.YLabel.String = {'v_{||}','(10^3 km/s)'};  
  set(hca,'ColorOrder',mms_colors('122'))
  %irf_legend(hca,{'v_{ph}'},[0.55 0.7],'fontsize',12);
  %irf_legend(hca,{'v_{trap}'},[0.55 0.99],'fontsize',12);
  %irf_legend(hca,{'v_{trap}'},[0.55 0.3],'fontsize',12);
  
  hsurf = findobj(hca.Children,'Type','Surface');
  hsurf.FaceAlpha = 1;
  hline = findobj(hca.Children,'Type','Line');
  c_eval('hline(?).LineWidth = 1.5;',1:numel(hline))
end

h(4).Position(4) = 2*h(4).Position(4);
hcb.Position(4) = h(4).Position(4)*0.95;
h(4).YLim = [-59.9 59.9];
h(4).YTick = [-40 -20 0 20 40];

%h(5).CLim = [-7 -3];
hca = irf_panel('phase velocity');
%hca.CLim = [-5 -2.5];
%h(1).CLim = [-5 -2.5];

%colormap(hca,'jet')
irf_zoom(h,'x',tint_fred)
irf_zoom(h(1:2),'y')
irf_plot_axis_align
%h(4).CLim = [-35 -28]+12
%colormap(cn.cmap('blue_white'));
colormap('jet')
delete(irf_panel('delete'))
