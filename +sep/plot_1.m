ic = 1;
event = 3;

%% Set datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/data/mms');
db_info = datastore('mms_db');   
localuser = datastore('local','user');

%% Tint
sep.get_tints;

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
c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',1:4);
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

%% Fred
eDist = ePDist1.tlim(tint_fred);
ve = gseVe1.tlim(eDist.time).resample(eDist);
 % keep in mind that this also affects the velocity at lower energies
scpot = scPot1.resample(eDist);
lowerelim = 100;
energies = sqrt(2*eDist.depend{1}(1,:)*units.eV/units.me)/1000; % km/s
vgmax = 70000;
%vg = [-energies(end:-1:1) 0 energies];
vg = -vgmax:1000:vgmax;
vg(abs(vg)>70000) = [];
nMC = 500;
tic; ef1D = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'scpot',scpot,'lowerelim',lowerelim,'vg',vg); toc % reduced distribution along B
%tic; ef1D_sep = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'scpot',scpot,'lowerelim',lowerelim,'vg',vg); toc % reduced distribution along B
%tic; ef1D_sheet = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'scpot',scpot,'lowerelim',lowerelim,'vg',vg); toc % reduced distribution along B
%tic; ef1D_lobe = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'scpot',scpot,'lowerelim',lowerelim,'vg',vg); toc % reduced distribution along B
disp('Done.')

%% Pick out parameters in the lobe, sheet, and separatrix intervals
c_eval('B_lobe = mean(gseB?.abs.tlim(tint_lobe).data,1);',ic)

c_eval('n_lobe = mean(ne?.tlim(tint_lobe).data,1);',ic)
c_eval('n_sheet = mean(ne?.tlim(tint_sheet).data,1);',ic)
c_eval('n_sep = mean(ne?.tlim(tint_sep).data,1);',ic)

c_eval('be_lobe = mean(beta?e.tlim(tint_lobe).data,1);',ic)
c_eval('be_sheet = mean(beta?e.tlim(tint_sheet).data,1);',ic)
c_eval('be_sep = mean(beta?e.tlim(tint_sep).data,1);',ic)

c_eval('Tepar_lobe = mean(facTe?.tlim(tint_lobe).data(:,1,1),1);',ic)
c_eval('Tepar_sheet = mean(facTe?.tlim(tint_sheet).data(:,1,1),1);',ic)
c_eval('Tepar_sep = mean(facTe?.tlim(tint_sep).data(:,1,1),1);',ic)

c_eval('Teperp_lobe  = 0.5*mean(facTe?.tlim(tint_lobe).data(:,2,2)+facTe?.tlim(tint_lobe).data(:,3,3),1);',ic)
c_eval('Teperp_sheet = 0.5*mean(facTe?.tlim(tint_sheet).data(:,2,2)+facTe?.tlim(tint_sheet).data(:,3,3),1);',ic)
c_eval('Teperp_sep   = 0.5*mean(facTe?.tlim(tint_sep).data(:,2,2)+facTe?.tlim(tint_sep).data(:,3,3),1);',ic)

wce0 = units.e*B_lobe*1e-9/units.me; % fce0*2pi
wpe0 = sqrt(n_sheet*1e6*units.e^2/units.eps0/units.me); % fpe0*2*pi

info_str = {sprintf('B_{lobe}=%.0f nT',B_lobe);...
            sprintf('n_{lobe}=%.3f cc',n_lobe);...
            sprintf('n_{sheet}=%.3f cc',n_sheet);...
            sprintf('n_{sep}=%.3f cc',n_sep);...
            sprintf('T_{e,par,lobe}=%.0f eV',Tepar_lobe);...
            sprintf('T_{e,par,sheet}=%.0f eV',Tepar_sheet);...
            sprintf('T_{e,par,sep}=%.0f eV',Tepar_sep);...
            sprintf('T_{e,perp,lobe}=%.0f eV',Teperp_lobe);...
            sprintf('T_{e,perp,sheet}=%.0f eV',Teperp_sheet);...
            sprintf('T_{e,perp,sep}=%.0f eV',Teperp_sep);...
            sprintf('f_{ce,0}=%.0f Hz',wce0/2/pi);...
            sprintf('f_{pe,0}=%.0f Hz',wpe0/2/pi);...
            sprintf('f_{pe,0}/f_{ce,0}=%.1f',wpe0/wce0);...
            sprintf('beta_{e,lobe}=%.3f',be_lobe);...
            sprintf('beta_{e,sheet}=%.3f',be_sheet);...
            sprintf('beta_{e,sep}=%.3f',be_sep);...
            };
colors = mms_colors('xzy');
info_color = [0 0 0; colors; colors; colors; 0 0 0; 0 0 0; 0 0 0; colors];

tint_day_utc = tint.utc; tint_day_utc = tint_day_utc(1,1:10);
tint_lobe_utc = tint_lobe.utc; tint_lobe_utc = tint_lobe_utc(:,12:23);
tint_sheet_utc = tint_sheet.utc; tint_sheet_utc = tint_sheet_utc(:,12:23);
tint_sep_utc = tint_sep.utc; tint_sep_utc = tint_sep_utc(:,12:23);
tint_phi_utc = tint_phi.utc; tint_phi_utc = tint_phi_utc(:,12:23);

str_print=sprintf('sep_ov_%s_%s_%s',tint_day_utc,tint_sep_utc(1,:),tint_sep_utc(2,:));
str_print(strfind(str_print,':'))=[];

%% Get wave parameters for tint_phi interval

%% Plot overview figure with focus on electrons
figure(201)
npanels_large = 10;
npanels_zoom = 1;
npanels = 12;
cmap = 'jet';
[h,h2] = initialize_combined_plot(npanels,3,2,0.4,'vertical'); % horizontal
ic = 1;
iisub = 0;
cmap = colormap('jet');
zoomy = [];
fig = gcf;
fig.Position(3) = fig.Position(4)*1.5;

fred_min = 1e-6;
fontsize = 10;

if 1 % e DEF omni 64
  iisub = iisub + 1;
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
if 1 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % neiisub = iisub + 1;
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 1 % Jx, moments and curl  iisub = iisub + 1;
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Jx fpi fgm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJcurl.x.tlim(tint+[1 0])},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J_x','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'FPI','FGM'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % J  iisub = iisub + 1;
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Vi  iisub = iisub + 1;
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve  
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve perp par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
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
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 0 % e DEF omni 64
  iisub = iisub + 1;
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
if 1 % Te par perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
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
if 1 % ePDist pa 64
  iisub = iisub + 1;
  hca = irf_panel('e pitch');  
  eint = [100 32000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec,''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1],'fontsize',fontsize)
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  hca.YLim = [0 180];
end
if 1 % e psd vpar
  iisub = iisub + 1;
  hca = irf_panel('fe reduced');
  fred_to_plot = ef1D; fred_to_plot.data(fred_to_plot.data < fred_min) = NaN;
  irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));  
  hca.YLim = fred_to_plot.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  %irf_legend(hca,[num2str(fred_to_plot.ancillary.vint(1),'%.0f') '<v_\perp<' num2str(fred_to_plot.ancillary.vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(lowerelim) ' eV'],[0.98 0.99],'color',0*[1 1 1],'fontsize',fontsize)
  irf_legend(hca,['f_{e} >' num2str(fred_min) ' s/m^4'],[0.98 0.70],'color',0*[1 1 1],'fontsize',fontsize)
end
if 0 % E
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x.tlim(tint),gseE?perp.y.tlim(tint),gseE?perp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
hca = irf_panel('delete');
if 1 % e psd vpar
  iisub = iisub + 1;
  hca = irf_panel('fe reduced zoomin');
  fred_to_plot = ef1D; fred_to_plot.data(fred_to_plot.data < fred_min) = NaN;
  irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));  
  hca.YLim = fred_to_plot.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'};   
  irf_legend(hca,['E_{e} >' num2str(lowerelim) ' eV'],[0.98 1.1],'color',0*[1 1 1],'fontsize',fontsize)
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels_large npanels_large + 1 + [1:npanels_zoom]]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end


irf_zoom(h(1:npanels_large),'x',tint)
irf_zoom(h(npanels_large+1+[1:npanels_zoom]),'x',tint_fred)
irf_zoom(h(zoomy),'y') 
irf_plot_axis_align
[hline1,hline2] = irf_plot_zoomin_lines_between_panels(h(npanels_large),h(npanels_large+2)); 
irf_timeaxis(h(npanels_large))

h(1).Title.String = sprintf('MMS %g, event %g',ic,event);

clear hmark_lobe hmark_sheet hmark_sep
c_eval('hmark_lobe(?) = irf_pl_mark(h(?),tint_lobe, mms_colors(''x''));',1:npanels_large)
c_eval('hmark_sheet(?) = irf_pl_mark(h(?),tint_sheet, mms_colors(''z''));',1:npanels_large)
c_eval('hmark_sep(?) = irf_pl_mark(h(?),tint_sep, mms_colors(''y''));',1:npanels_large)

hca = h(2);
set(hca,'ColorOrder',info_color)
irf_legend(hca,info_str,[1.01 0.99],'fontsize',fontsize)
for ii = 1:npanels
  h(ii).FontSize = 12;
end

hca = irf_panel('e DEF omni 64'); hca.XGrid = 'off'; hca.YGrid = 'off';
delete(irf_panel('delete'))

%% Plot acceleration potential and electron particle distributions
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?.resample(gseVe?);',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
%c_eval('ePitch = ePitch?;',ic)

times_phi_utc = ef1D.tlim(tint_phi).time.utc; times_phi_utc = times_phi_utc(:,12:23);

fontsize = 8;
  
colors = mms_colors('xzy');
isub = 1;
if 1 % average fred during sep, lobe and sheet interval
  hca = h2(isub); isub = isub + 1;
  
  fered_sep = ef1D.tlim(tint_sep);
  plot_v_sep = nanmean(fered_sep.depend{1},1);
  plot_f_sep = nanmean(fered_sep.data,1);
  
  fered_lobe = ef1D.tlim(tint_lobe);
  plot_v_lobe = nanmean(fered_lobe.depend{1},1);
  plot_f_lobe = nanmean(fered_lobe.data,1);
  
  fered_sheet = ef1D.tlim(tint_sheet);
  plot_v_sheet = nanmean(fered_sheet.depend{1},1);
  plot_f_sheet = nanmean(fered_sheet.data,1);
  
  h_fered = plot(hca,plot_v_lobe*1e-3,plot_f_lobe,...
                     plot_v_sheet*1e-3,plot_f_sheet,...
                     plot_v_sep*1e-3,plot_f_sep);
  c_eval('h_fered(?).Color = colors(?,:);',1:3)
  hca.XLabel.String = 'v (10^3 km/s)';
  hca.YLabel.String = sprintf('f_e (%s)',fered_sep.units);
  hca.XLim = [-40 40];
  
  set(hca,'colororder',colors);
  tint_str = {sprintf('%s - %s',tint_lobe_utc(1,:),tint_lobe_utc(2,:));...
              sprintf('%s - %s',tint_sheet_utc(1,:),tint_sheet_utc(2,:));...
              sprintf('%s - %s',tint_sep_utc(1,:),tint_sep_utc(2,:))};
  irf_legend(hca,tint_str,[0.01 0.98],'fontsize',fontsize)
  hca.Title.String = sprintf('E > %g eV',lowerelim);
end

%ndist_to_combine = 3;
n_int = 3;
lowerelim_peaks = 400;
Elim_axis = 4000;
colors = mms_colors('matlab');
clear plot_f_phi f_sep_utc hmark_all str_Emax Emax Emax_within_interval 
clear Emax_av Emax_av_ind plot_f_phi_Emax Emax_utc f_phi_utc

if 1 % divide f into intervals, common for both v,E x-axis, mark TSeries
  fred_phi = ef1D.tlim(tint_phi);
  plot_v_phi = mean(fred_phi.depend{1},1);
  plot_E_phi = sign(plot_v_phi)*0.5*units.me.*(plot_v_phi*1e3).^2/units.e;   
  
  fred_phi_peaks = fred_phi;   
  fred_phi_peaks.data(abs(repmat(plot_E_phi,fred_phi_peaks.length,1))<lowerelim_peaks) = 0;
  [val,ind] = max(fred_phi_peaks.data');
  Emax = (plot_E_phi(ind));
     
  nf = fred_phi.length;
  [ind1,ind2] = ind12(nf,n_int);
  for i_int = 1:n_int
    f_phi_utc{i_int,1} = sprintf('%s - %s',times_phi_utc(ind1(i_int),:),times_phi_utc(ind2(i_int),:));
    Emax_tmp = Emax(ind1(i_int):ind2(i_int));
    [val_tmp_Emax,ind_tmp_Emax] = max(abs(Emax_tmp));
    Emax_within_interval(i_int) = Emax_tmp(ind_tmp_Emax);
    %Emax_utc = 
        
    plot_f_phi(i_int,:) = mean(fred_phi.data(ind1(i_int):ind2(i_int),:),1);    
    plot_f_phi_Emax(i_int,:) = fred_phi.data(ind1(i_int)+ind_tmp_Emax-1,:);
    
    Emax_tmp = plot_f_phi(i_int,:); Emax_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(Emax_tmp);
    Emax_av(i_int) = (plot_E_phi(ind));
    Emax_av_ind(i_int) = ind;
    f_at_Emax_av(i_int) = plot_f_phi(i_int,ind);
    
    str_Emax{i_int,1} = sprintf('phimax =[%.0f, %.0f] eV',Emax_av(i_int),Emax_within_interval(i_int));
    
    [hmark,cc] = irf_pl_mark(h(npanels_large + 2),fred_phi.time([ind1(i_int) ind2(i_int)]).epochUnix');    
    if isa(hmark,'matlab.graphics.primitive.Line')
      hmark.Color = colors(i_int,:);
    elseif isa(hmark,'matlab.graphics.primitive.Patch')
      hmark.YData = hmark.YData + 120;
      hmark.ZData = [1 1 1 1];
      hmark.EdgeColor = [0 0 0];
      hmark.FaceColor = colors(i_int,:);
    end    
    hmark_all{i_int} = hmark;
  end  
end
if 1 % all fred during sep interval, vs vpar
  hca = h2(isub); isub = isub + 1;
  h_fered = plot(hca,plot_v_sep*1e-3,plot_f_phi);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)  
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = sprintf('f_e (%s)',fered_sep.units);
  hca.XLim = [-40 40];
  set(hca,'colororder',colors)  
  hca.Title.String = 'Average over indicated tint';
  if plot_E_phi(ind) > 0
    irf_legend(hca,f_phi_utc,[0.01 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,f_phi_utc,[0.98 0.98],'fontsize',fontsize)
  end
end
if 1 % all fred during sep interval, vs +-Epar
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi*1e0,plot_f_phi);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    h_fered_star = plot(hca,Emax_av(i_int),f_at_Emax_av(i_int),'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.XLim = Elim_axis*[-1 1]+ 0.5*mean(Emax_av);
  hca.XTick = -30000:1000:30000;
     
  hca.Title.String = {'Average over indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};

  if 1 % plot lowerelim_peaks
    hold(hca,'on')
    plot(hca,lowerelim_peaks*[-1 -1],hca.YLim,'--k',lowerelim_peaks*[1 1],hca.YLim,'--k')
    hold(hca,'off')
  end
  
  set(hca,'colororder',colors)
  if plot_E_phi(ind) > 0
    irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
  end
end
if 1 % all fred during sep interval, vs +-Epar, max within interal
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi,plot_f_phi_Emax);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    f_tmp = plot_f_phi_Emax(i_int,:);
    f_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(abs(f_tmp));
    h_fered_star = plot(hca,plot_E_phi(ind),val,'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.XLim = Elim_axis*[-1 1]+ 0.5*mean(Emax);
  hca.XTick = -30000:1000:30000;
  if hca.XLim(2)<max(Emax)
    hca.XLim = max(Emax)*1*[-1 1] + 0.5*max(Emax);
    hca.XTick = -30000:2000:30000;
  end
     
  hca.Title.String = sprintf('E > %g eV',lowerelim_peaks);
  if 1 % plot lowerelim_peaks
    hold(hca,'on')
    plot(hca,lowerelim_peaks*[-1 -1],hca.YLim,'--k',lowerelim_peaks*[1 1],hca.YLim,'--k')
    hold(hca,'off')
  end
  
  set(hca,'colororder',colors)
  if plot_E_phi(ind) > 0
    irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
  end
  hca.Title.String = {'Max phi in indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};
end

n_intervals_sep = 5;
vgmax = 60000;
vg = -vgmax:2000:vgmax;
%lowerelim = 40;
if 1 % 1 reduced distribution, 1st half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[0 1];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  %axis(hca,'equal')
  axis(hca,'square')
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
end
if 1 % 1 reduced distribution, 2nd half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[1 2];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  

  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  %axis(hca,'equal')
  axis(hca,'square')
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
end
if 0 % 1 reduced distribution, 3rd third of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[2 3];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
  %axis(hca,'equal')
  axis(hca,'square')
end

for ih = 1:numel(h2)
  h2(ih).FontSize = fontsize;
end
%for ih = 4:6, h2(ih).CLim = [-15 -10.2]; end
  
%% Plot overview figure with focus on waves
figure(202)
npanels_large = 5;
npanels_zoom = 2;
npanels_delete = 1;
npanels = npanels_large + npanels_zoom + npanels_delete;
cmap = 'jet';
[h,h2] = initialize_combined_plot(npanels,3,2,0.4,'vertical'); % horizontal
ic = 1;
iisub = 0;
cmap = colormap('jet');
zoomy = [];
fig = gcf;
fig.Position(3) = fig.Position(4)*1.5;

fred_min = 1e-6;
fontsize = 10;

if 1 % e DEF omni 64
  iisub = iisub + 1;
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
if 0 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % neiisub = iisub + 1;
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 0 % Jx, moments and curl  iisub = iisub + 1;
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Jx fpi fgm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJcurl.x.tlim(tint+[1 0])},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J_x','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'FPI','FGM'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % J  iisub = iisub + 1;
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi  iisub = iisub + 1;
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve  
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve perp par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
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
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 0 % e DEF omni 64
  iisub = iisub + 1;
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
if 1 % Te par perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
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
if 0 % ePDist pa 64
  iisub = iisub + 1;
  hca = irf_panel('e pitch');  
  eint = [100 32000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec,''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1],'fontsize',fontsize)
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  hca.YLim = [0 180];
end
if 0 % e psd vpar
  iisub = iisub + 1;
  hca = irf_panel('fe reduced');
  fred_to_plot = ef1D; fred_to_plot.data(fred_to_plot.data < fred_min) = NaN;
  irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));  
  hca.YLim = fred_to_plot.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  %irf_legend(hca,[num2str(fred_to_plot.ancillary.vint(1),'%.0f') '<v_\perp<' num2str(fred_to_plot.ancillary.vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(lowerelim) ' eV'],[0.98 0.99],'color',0*[1 1 1],'fontsize',fontsize)
  irf_legend(hca,['f_{e} >' num2str(fred_min) ' s/m^4'],[0.98 0.70],'color',0*[1 1 1],'fontsize',fontsize)
end
if 0 % E
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x.tlim(tint),gseE?perp.y.tlim(tint),gseE?perp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
for ii = 1:npanels_delete
  hca = irf_panel(sprintf('delete %g',ii));
  iisub = iisub + 1;
end
if 1 % e psd vpar
  iisub = iisub + 1;
  hca = irf_panel('fe reduced zoomin');
  fred_to_plot = ef1D; fred_to_plot.data(fred_to_plot.data < fred_min) = NaN;
  irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));  
  hca.YLim = fred_to_plot.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'};   
  irf_legend(hca,['E_{e} >' num2str(lowerelim) ' eV'],[0.98 1.1],'color',0*[1 1 1],'fontsize',fontsize)
end
if 1 % E par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par zoomin');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels_large npanels-npanels_delete-1+[1:npanels_zoom]]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end


irf_zoom(h(1:npanels_large),'x',tint)
irf_zoom(h(npanels_large + npanels_delete + [1:npanels_zoom]),'x',tint_fred)
irf_zoom(h(zoomy),'y') 
irf_plot_axis_align
[hline1,hline2] = irf_plot_zoomin_lines_between_panels(h(npanels_large),h(npanels_large + npanels_delete + 1)); 
irf_timeaxis(h(npanels_large))

h(1).Title.String = sprintf('MMS %g, event %g',ic,event);

clear hmark_lobe hmark_sheet hmark_sep
c_eval('hmark_lobe(?) = irf_pl_mark(h(?),tint_lobe, mms_colors(''x''));',1:npanels_large)
c_eval('hmark_sheet(?) = irf_pl_mark(h(?),tint_sheet, mms_colors(''z''));',1:npanels_large)
c_eval('hmark_sep(?) = irf_pl_mark(h(?),tint_sep, mms_colors(''y''));',1:npanels_large)

hca = h(2);
set(hca,'ColorOrder',info_color)
irf_legend(hca,info_str,[1.01 0.99],'fontsize',fontsize)
for ii = 1:npanels
  h(ii).FontSize = 12;
end

hca = irf_panel('e DEF omni 64'); hca.XGrid = 'off'; hca.YGrid = 'off';

for ii = 1:npanels_delete
  delete(irf_panel(sprintf('delete %g',ii)));  
end

%% Plot electron particle distributions
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?.resample(gseVe?);',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
%c_eval('ePitch = ePitch?;',ic)

times_phi_utc = ef1D.tlim(tint_phi).time.utc; times_phi_utc = times_phi_utc(:,12:23);

fontsize = 8;
  
colors = mms_colors('xzy');
isub = 1;
if 1 % average fred during sep, lobe and sheet interval
  hca = h2(isub); isub = isub + 1;
  
  fered_sep = ef1D.tlim(tint_sep);
  plot_v_sep = nanmean(fered_sep.depend{1},1);
  plot_f_sep = nanmean(fered_sep.data,1);
  
  fered_lobe = ef1D.tlim(tint_lobe);
  plot_v_lobe = nanmean(fered_lobe.depend{1},1);
  plot_f_lobe = nanmean(fered_lobe.data,1);
  
  fered_sheet = ef1D.tlim(tint_sheet);
  plot_v_sheet = nanmean(fered_sheet.depend{1},1);
  plot_f_sheet = nanmean(fered_sheet.data,1);
  
  h_fered = plot(hca,plot_v_lobe*1e-3,plot_f_lobe,...
                     plot_v_sheet*1e-3,plot_f_sheet,...
                     plot_v_sep*1e-3,plot_f_sep);
  c_eval('h_fered(?).Color = colors(?,:);',1:3)
  hca.XLabel.String = 'v (10^3 km/s)';
  hca.YLabel.String = sprintf('f_e (%s)',fered_sep.units);
  hca.XLim = [-40 40];
  
  set(hca,'colororder',colors);
  tint_str = {sprintf('%s - %s',tint_lobe_utc(1,:),tint_lobe_utc(2,:));...
              sprintf('%s - %s',tint_sheet_utc(1,:),tint_sheet_utc(2,:));...
              sprintf('%s - %s',tint_sep_utc(1,:),tint_sep_utc(2,:))};
  irf_legend(hca,tint_str,[0.01 0.98],'fontsize',fontsize)
  hca.Title.String = sprintf('E > %g eV',lowerelim);
end

%ndist_to_combine = 3;
n_int = 3;
lowerelim_peaks = 400;
Elim_axis = 4000;
colors = mms_colors('matlab');
clear plot_f_phi f_sep_utc hmark_all str_Emax Emax Emax_within_interval 
clear Emax_av Emax_av_ind plot_f_phi_Emax Emax_utc f_phi_utc

if 1 % divide f into intervals, common for both v,E x-axis, mark TSeries
  fred_phi = ef1D.tlim(tint_phi);
  plot_v_phi = mean(fred_phi.depend{1},1);
  plot_E_phi = sign(plot_v_phi)*0.5*units.me.*(plot_v_phi*1e3).^2/units.e;   
  
  fred_phi_peaks = fred_phi;   
  fred_phi_peaks.data(abs(repmat(plot_E_phi,fred_phi_peaks.length,1))<lowerelim_peaks) = 0;
  [val,ind] = max(fred_phi_peaks.data');
  Emax = (plot_E_phi(ind));
     
  nf = fred_phi.length;
  [ind1,ind2] = ind12(nf,n_int);
  for i_int = 1:n_int
    f_phi_utc{i_int,1} = sprintf('%s - %s',times_phi_utc(ind1(i_int),:),times_phi_utc(ind2(i_int),:));
    Emax_tmp = Emax(ind1(i_int):ind2(i_int));
    [val_tmp_Emax,ind_tmp_Emax] = max(abs(Emax_tmp));
    Emax_within_interval(i_int) = Emax_tmp(ind_tmp_Emax);
    %Emax_utc = 
        
    plot_f_phi(i_int,:) = mean(fred_phi.data(ind1(i_int):ind2(i_int),:),1);    
    plot_f_phi_Emax(i_int,:) = fred_phi.data(ind1(i_int)+ind_tmp_Emax-1,:);
    
    Emax_tmp = plot_f_phi(i_int,:); Emax_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(Emax_tmp);
    Emax_av(i_int) = (plot_E_phi(ind));
    Emax_av_ind(i_int) = ind;
    f_at_Emax_av(i_int) = plot_f_phi(i_int,ind);
    
    str_Emax{i_int,1} = sprintf('phimax =[%.0f, %.0f] eV',Emax_av(i_int),Emax_within_interval(i_int));
    
    [hmark,cc] = irf_pl_mark(h(npanels_large + 2),fred_phi.time([ind1(i_int) ind2(i_int)]).epochUnix');    
    if isa(hmark,'matlab.graphics.primitive.Line')
      hmark.Color = colors(i_int,:);
    elseif isa(hmark,'matlab.graphics.primitive.Patch')
      hmark.YData = hmark.YData + 120;
      hmark.ZData = [1 1 1 1];
      hmark.EdgeColor = [0 0 0];
      hmark.FaceColor = colors(i_int,:);
    end    
    hmark_all{i_int} = hmark;
  end  
end
if 1 % all fred during sep interval, vs vpar
  hca = h2(isub); isub = isub + 1;
  h_fered = plot(hca,plot_v_sep*1e-3,plot_f_phi);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)  
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = sprintf('f_e (%s)',fered_sep.units);
  hca.XLim = [-40 40];
  set(hca,'colororder',colors)  
  hca.Title.String = 'Average over indicated tint';
  if plot_E_phi(ind) > 0
    irf_legend(hca,f_phi_utc,[0.01 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,f_phi_utc,[0.98 0.98],'fontsize',fontsize)
  end
end
if 1 % all fred during sep interval, vs +-Epar
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi*1e0,plot_f_phi);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    h_fered_star = plot(hca,Emax_av(i_int),f_at_Emax_av(i_int),'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.XLim = Elim_axis*[-1 1]+ 0.5*mean(Emax_av);
  hca.XTick = -30000:1000:30000;
     
  hca.Title.String = {'Average over indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};

  if 1 % plot lowerelim_peaks
    hold(hca,'on')
    plot(hca,lowerelim_peaks*[-1 -1],hca.YLim,'--k',lowerelim_peaks*[1 1],hca.YLim,'--k')
    hold(hca,'off')
  end
  
  set(hca,'colororder',colors)
  if plot_E_phi(ind) > 0
    irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
  end
end
if 1 % all fred during sep interval, vs +-Epar, max within interal
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi,plot_f_phi_Emax);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    f_tmp = plot_f_phi_Emax(i_int,:);
    f_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(abs(f_tmp));
    h_fered_star = plot(hca,plot_E_phi(ind),val,'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.XLim = Elim_axis*[-1 1]+ 0.5*mean(Emax);
  hca.XTick = -30000:1000:30000;
  if hca.XLim(2)<max(Emax)
    hca.XLim = max(Emax)*1*[-1 1] + 0.5*max(Emax);
    hca.XTick = -30000:2000:30000;
  end
     
  hca.Title.String = sprintf('E > %g eV',lowerelim_peaks);
  if 1 % plot lowerelim_peaks
    hold(hca,'on')
    plot(hca,lowerelim_peaks*[-1 -1],hca.YLim,'--k',lowerelim_peaks*[1 1],hca.YLim,'--k')
    hold(hca,'off')
  end
  
  set(hca,'colororder',colors)
  if plot_E_phi(ind) > 0
    irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
  end
  hca.Title.String = {'Max phi in indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};
end

n_intervals_sep = 2;
vgmax = 60000;
vg = -vgmax:2000:vgmax;
%lowerelim = 40;
if 1 % 1 reduced distribution, 1st half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[0 1];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  %axis(hca,'equal')
  axis(hca,'square')
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
end
if 1 % 1 reduced distribution, 2nd half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[1 2];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  %axis(hca,'equal')
  axis(hca,'square')
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
end
if 0 % 1 reduced distribution, 3rd third of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[2 3];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
  %axis(hca,'equal')
  axis(hca,'square')
end

for ih = 1:numel(h2)
  h2(ih).FontSize = fontsize;
end
%for ih = 4:6, h2(ih).CLim = [-15 -10.2]; end
  