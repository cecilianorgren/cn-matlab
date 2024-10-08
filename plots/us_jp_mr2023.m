% US-Japan
%% Load data
ic = [3];
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854

% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db','/Users/cno062/Data/MMS');
db_info = datastore('mms_db');

units = irf_units;

% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);

%c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)

% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
%c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',tint); toc',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);

% Load spacecraft position
disp('Loading spacecraft position...')
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4);

% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);

% Particle moments
% Skymap distributions
if 1
  %%
disp('Loading skymaps...')
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
if 0
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic)
c_eval('iPDist?_onecount = iPDist?; iPDist?_onecount.data = (iPDist?_onecount.data./iPDistErr?.data).^2;',ic)
end
end
% Pressure and temperature
disp('Loading pressure and temperature...'); tic
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic); toc

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
disp('Loading density...'); tic;
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic); toc

% Velocity
disp('Loading bulk velocities...'); tic
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); toc
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic); toc

% Calculate stress tensor from pressure, density and speed
% P_psd(:,1,1) = P_psd(:,1,1)-pmass*n_psd.*V_psd(:,1).*V_psd(:,1);
comps = ['x','y','z'];
c_eval('Se? = gsePe?.data*0;',ic);
c_eval('Si? = gsePi?.data*0;',ic);
for ic1 = 1:3
  for ic2 = 1:3
    
    c_eval('Se?(:,ic1,ic2) = gsePe?.data(:,ic1,ic2)*1e-9 + units.me*1e6*1e3*1e3*ne?.data.*gseVe?.data(:,ic1).*gseVe?.data(:,ic2);',ic);
    c_eval('Si?(:,ic1,ic2) = gsePi?.data(:,ic1,ic2)*1e-9 + units.mp*1e6*1e3*1e3*ni?.data.*gseVi?.data(:,ic1).*gseVi?.data(:,ic2);',ic);
    c_eval('gseSe? = irf.ts_tensor_xyz(gsePe?.time,Se?*1e9); gseSe?.name = ''Stress tensor''; gseSe?.units = ''nPa'';',ic)
    c_eval('gseSi? = irf.ts_tensor_xyz(gsePi?.time,Si?*1e9); gseSi?.name = ''Stress tensor''; gseSi?.units = ''nPa'';',ic)
  end
end
%




c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s

% EDR signatures
c_eval('facPepp? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('facPeqq? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic);
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic);

% Compute Q and Dng from facPepp
c_eval('Q? = (facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2)./(facPepp?.yy.data.^2+2*facPepp?.yy.data.*facPepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
c_eval('Dng? = sqrt(8*(facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2))./(facPepp?.xx.data+2*facPepp?.yy.data);',ic);
c_eval('Dng? = irf.ts_scalar(ne?.time,Dng?);',ic);

% Compute agyrotropy Aphi from facPeqq
c_eval('agyro? = 2*(facPeqq?.yy-facPeqq?.zz)/(facPeqq?.yy+facPeqq?.zz); agyro? = agyro?.abs',ic);

% Compute temperature ratio An
c_eval('Temprat? = facPepp?.xx/(facPepp?.yy);',ic);

c_eval('Te?par = facPepp?.xx;',ic);
c_eval('Te?perp = facPepp?.yy;',ic);

% Rotate things into new coordinate system 
L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9230];
lmn = [L;M;N];


% Rotate data
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)

c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaVe? = gseVe?*lmn''; mvaVe?.name = ''Ve LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaVExB? = gseVExB?*lmn''; mvaVExB?.name = ''VExB LMN'';',ic)
%c_eval('mvaJ? = gseJ?*lmn''; mvaJ?.name = ''J LMN'';',ic)
%c_eval('mvaJe? = gseJe?*lmn''; mvaJe?.name = ''Je LMN'';',ic)
%c_eval('mvaJi? = gseJi?*lmn''; mvaJi?.name = ''Ji LMN'';',ic)
%mvaJcurl = gseJcurl*lmn'; mvaJcurl.name = 'J LMN CURL';
c_eval('mvaPi? = lmn*gsePi?*lmn''; mvaPi?.units = gsePi?.units;',ic)
c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;',ic)
c_eval('mvaTi? = lmn*gseTi?*lmn''; mvaTi?.units = gseTi?.units;',ic)
c_eval('mvaTe? = lmn*gseTe?*lmn''; mvaTe?.units = gseTe?.units;',ic)


c_eval('mvaSi? = lmn*gseSi?*lmn''; mvaSi?.units = gseSi?.units;',ic)
c_eval('mvaSe? = lmn*gseSe?*lmn''; mvaSe?.units = gseSe?.units;',ic)

c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;',ic)


c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';',1:4)

%% Spacecraft constellation
tt = irf_time('2017-07-11T22:33:58.000000000Z','utc>EpochTT');
R0 = (mvaR1 + mvaR2.resample(mvaR1) + mvaR3.resample(mvaR1) + mvaR4.resample(mvaR1))/4;
c_eval('r? = mvaR?-R0;',1:4)
c_eval('r? = r?.resample(tt).data;',1:4)
%c_eval('dr(?,!) = dr? - dr!;',1:4,1:4)
r = [r1;r2;r3;r4];
rmax = max(r);
rmin = min(r);

h = subplot(1,1,1);
isub = 1;
hca = h(isub);
symbols = {'o','s','^','<'};
colors = mms_colors('1234');
colors(1,:) = [0.2 0.2 0.2];

colors = [0.2000    0.2000    0.2000;
          0.9000    0.2000         0;
               0    0.8000         0;
          0.1000    0.4000    0.9000];
linewidth = 2;
markersize = 20;

holdon = 0;
for isc = 1:4
  plot3(hca,r(isc,1),r(isc,2),r(isc,3),'s','markersize',markersize,'linewidth',linewidth,'color',colors(isc,:));
  if not(holdon)
    hold(hca,'on')
  end
end

hca.XLim = [-15 15];
hca.YLim = [-15 15];
hca.ZLim = [-15 15];

hca.XLim = [rmin(1) rmax(1)]+[-1 1]*2;
hca.YLim = [rmin(2) rmax(2)]+[-1 1]*2;
hca.ZLim = [rmin(3) rmax(3)]+[-1 1]*2;

xmin = hca.XLim(2);
ymin = hca.YLim(2);
zmin = hca.ZLim(1);
for isc = 1:4
  plot3(hca,xmin*[1 1],r(isc,2),r(isc,3),'s','markersize',10,'linewidth',1,'color',colors(isc,:).^0.125);  
  plot3(hca,r(isc,1),ymin*[1 1],r(isc,3),'s','markersize',10,'linewidth',1,'color',colors(isc,:).^0.125);  
  plot3(hca,r(isc,1),r(isc,2),zmin*[1 1],'s','markersize',10,'linewidth',1,'color',colors(isc,:).^0.125);  
end

for i1 = 1:4
  for i2 = i1:4
    plot3(hca,[r(i1,1) r(i2,1)],[r(i1,2) r(i2,2)],[r(i1,3) r(i2,3)],'linewidth',1,'color',[0 0 0]);
  end
end
for i1 = 1:4
  for i2 = i1:4
    plot3(hca,[xmin xmin],[r(i1,2) r(i2,2)],[r(i1,3) r(i2,3)],':','linewidth',1,'color',[1 1 1]*0.7);
    plot3(hca,[r(i1,1) r(i2,1)],[ymin ymin],[r(i1,3) r(i2,3)],':','linewidth',1,'color',[1 1 1]*0.7);
    plot3(hca,[r(i1,1) r(i2,1)],[r(i1,2) r(i2,2)],[zmin zmin],':','linewidth',1,'color',[1 1 1]*0.7);
  end
end

hold(hca,'off')

hca.XLabel.String = 'L (km)';
hca.YLabel.String = 'M (km)';
hca.ZLabel.String = 'N (km)';

c_eval('h(?).FontSize = 14;',1:numel(h))

hca.Box = 'on';
%axis(hca,'square');

%daspect(hca,[1 1 1])
hca.DataAspectRatio = [1 1 1];
%hca.XLim = [rmin(1) rmax(1)]+[-1 1]*2;
%hca.YLim = [rmin(2) rmax(2)]+[-1 1]*2;
%hca.ZLim = [rmin(3) rmax(3)]+[-1 1]*2;


%% Figure: Overview 1
ic = 3;

tint_edr = irf.tint('2017-07-11T22:33:58.00Z/2017-07-11T22:34:10.00Z'); %20151112071854

npanels = 6;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
fontsize = 16;

isub = 0;
zoomy = [];

if 1 % B LMN
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize);
end 
if 1 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize);
end
if 1 % Pe-off LMN
  hca = irf_panel('Pe diag LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPe?.xx.tlim(tint)*1e3,mvaPe?.yy.tlim(tint)*1e3,mvaPe?.zz.tlim(tint)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'P_e','(pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'LL','MM','NN'}',[1.02 0.9],'fontsize',fontsize);
end
if 1 % Pe-off LMN
  hca = irf_panel('Pe off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPe?.xy.tlim(tint)*1e3,mvaPe?.xz.tlim(tint)*1e3,mvaPe?.yz.tlim(tint)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'P_e','(pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
end
if 1 % Pe-off LMN
  hca = irf_panel('Pe off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPe?.xy.tlim(tint)*1e3,mvaPe?.xz.tlim(tint)*1e3,mvaPe?.yz.tlim(tint)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'P_e','(pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
end
if 1 % Tpar, Tperp
  hca = irf_panel('Tepar, Teperp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Te?par,Te?perp},''comp'');',ic)  
  hca.YLabel.String = {'T_e','eV'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[1.02 0.9],'fontsize',fontsize);
end
if 1 % Non-gyro
  hca = irf_panel('Non-gyrotropies');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{agyro?,Dng?,Q?.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'','Non-gyrotropy'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'A_{\phi}','Dng','Q^{1/2}'}',[1.02 0.9],'fontsize',fontsize);
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.98],'color',[0 0 0],'fontsize',fontsize)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint_edr)
%irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))
h(end).XTickLabelRotation = 0;


%% 2D distribution and pressure contributions, MN, X times, with locations shown
ic = 3;
%t1 =  irf_time('2017-07-11T22:34:01.30Z','utc>EpochTT');
%t2 =  irf_time('2017-07-11T22:34:03.30Z','utc>EpochTT');
times = irf_time(['2017-07-11T22:34:01.30Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.03;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.4;
times = times(1):0.5:times(2);

times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z']);
%times = times([2 3 4])+0.25;
%times = times + 0.25;
times = times + 0.30;
%times = times + 0.06;
dt_dist = 2*0.062; % for two distributions

vint = [-Inf Inf];
vint = [-20000 20000];
vg = -60000:2000:60000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(3,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(1); 
h1.Position = position; 
%delete(h1_pos);
clear h2
for ip = 1:(nt*2)
  h2(ip) = subplot(3,nt,nt+ip);
end

hca = h1(1);
c_eval('mvaPe = mvaPe?;',ic)  
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaPe.xy*1e3,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaPe1.xy*1e3,mvaB1.abs},'comp')
hca.YLabel.String = 'P_{eLM} (pPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
h1.Position(4) = 0.16;
ppsum = [];
hp.EdgeColor = [0 0 0];
hp.FaceColor = [0.5 0.5 0.5];

times_exact = {};
isub = 1;
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca);
  hc_.Colorbar.YLabel.String = 'log_{10} f_e (s^2/m^5)';
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';

  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_distx = 0.030;
    B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(2)/b(1);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
    %B_
    irf_legend(hca,sprintf('B = %.1f nT', sqrt(sum(B_.^2))),[0.02 0.98],'color','k','fontsize',14)
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.x.resample(dist),ve.y.resample(dist));  
  hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';
  
  if 1
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    pptmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    ppsum(itime) = pptmp;
    
    if pptmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
  end

  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_dist = 0.030;
    B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(2)/b(1);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end
end

c_eval('hmark(?) = irf_pl_mark(h1(1),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]);',1:times.length)

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

h2(1).XLim = 0.99*vg([1 end])*1e-3;
h2(1).YLim = 0.99*vg([1 end])*1e-3;

c_eval('h2(?).XTick = -60:20:60; h2(?).YTick = -60:20:60;',1:numel(h2))

c_eval('h2(?).LineWidth = 1;',1:numel(h2))
c_eval('h2(?).FontSize = 14;',1:numel(h2))
c_eval('axis(h2(?),''square'');',1:numel(h2))

c_eval('h1(?).LineWidth = 1;',1:numel(h1))
c_eval('h1(?).FontSize = 14;',1:numel(h1))

%compact_panels(h2,0.0,00)
compact_panels(h2,0.005,00.005)
hb = findobj(gcf,'type','colorbar'); 
c_eval('hb(?).FontSize = 14;',1:numel(h2))
hb = hb(end:-1:1);
%
ihsub = [1:nt-1 nt+1:(2*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)

h1(1).Title.String = sprintf('MMS %g',ic);
h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
h1.Position(3) = h1.Position(3)*0.8;


if 1
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  tsPP = irf.ts_scalar(times,ppsum);
  hold(h1(1),'on')
  hpp = irf_plot(h1(1),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(1),'off')
end

%% 2D distribution and pressure contributions, LM, X times, comparing different satellites
ics = [3 4];
%t1 =  irf_time('2017-07-11T22:34:01.30Z','utc>EpochTT');
%t2 =  irf_time('2017-07-11T22:34:03.30Z','utc>EpochTT');
times = irf_time(['2017-07-11T22:34:01.30Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.03;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.4;
times = times(1):0.5:times(2);

times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z';...
  '2017-07-11T22:34:03.670000000Z';...
  '2017-07-11T22:34:04.300000000Z']);
%times = times([2 3 4])+0.25;
times = times + 0.05 + 0*0.25;
times = times + 0.05 + 2*0.03;


vint = [-Inf Inf];
vint = [-20000 20000];
vg = -60000:2000:60000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(5,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(1); 
h1.Position = position; 

%delete(h1_pos);
for ip = 1:(nt*2*2)
  h2(ip) = subplot(5,nt,nt+ip);
end

hca = h1(1);
set(hca,'colororder',mms_colors('12'))
c_eval('mvape_1 = mvaPe?;',ics(1))  
c_eval('mvape_2 = mvaPe?;',ics(2))  
irf_plot(hca,{mvape_1.yz,mvape_2.yz},'comp')
hold(hca,'on')
fhigh = 2;
irf_plot(hca,{mvape_1.yz.filt(0,fhigh,[],5),mvape_2.yz.filt(0,fhigh,[],5)},'comp')
hold(hca,'off')

hca.YLabel.String = 'P_{eMN} (nPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
hmark = irf_pl_mark(hca,times.epochUnix,'k');
c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
h1.Position(4) = 0.16;
%
isub = 1;
for ic = ics
  for itime = 1:times.length
    %hca = h2(isub); isub = isub + 1;
    time = times(itime);
    % Reduce distributions
    c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
    c_eval('scpot = scPot?.resample(dist);',ic)  
    vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
      
    % Plot
    hca = h2(isub); isub = isub + 1;
    [ha_,hb_,hc_] = vdf.plot_plane(hca);
    hc_.Colorbar.YLabel.String = 'log_{10} f_e(v_M,v_N) (s^2/m^5)';
    colormap(hca,pic_colors('candy4'))   
    hca.XLabel.String = 'v_M (10^3 km/s)';
    hca.YLabel.String = 'v_N (10^3 km/s)';
  end
  for itime = 1:times.length
    %hca = h2(isub); isub = isub + 1;
    time = times(itime);
    % Reduce distributions
    c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
    c_eval('scpot = scPot?.resample(dist);',ic)  
    c_eval('ve = mvaVe?;',ic)  
    vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
      
    % Plot
    hca = h2(isub); isub = isub + 1;
    [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.y.resample(dist),ve.z.resample(dist));  
    hc_.Colorbar.YLabel.String = 'f_e(v_M,v_N)(v_M-v_M^{bulk})(v_N-v_N^{bulk}) (1/m^3)';
    colormap(hca,pic_colors('blue_red'))
    hca.CLim = max(abs(hca.CLim))*[-1 1];
    hca.XLabel.String = 'v_M (10^3 km/s)';
    hca.YLabel.String = 'v_N (10^3 km/s)';
  end
end

hlinks_all = linkprop(h2,{'XLim','YLim'});
%hlinks_f = linkprop(h2([1:nt]*[0 2]'),{'CLim'});
%hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});
hlinks_f = linkprop(h2([1:nt 2*nt+(1:nt)]),{'CLim'});
hlinks_p = linkprop(h2([(nt+1):2*nt 2*nt+((nt+1):2*nt)]),{'CLim'});

h2(1).XLim = 0.99*vg([1 end])*1e-3;
h2(1).YLim = 0.99*vg([1 end])*1e-3;

c_eval('h2(?).XTick = -60:20:60; h2(?).YTick = -60:20:60;',1:numel(h2))

c_eval('h2(?).LineWidth = 1;',1:numel(h2))
c_eval('h2(?).FontSize = 14;',1:numel(h2))
c_eval('axis(h2(?),''square'');',1:numel(h2))

c_eval('h1(?).LineWidth = 1;',1:numel(h1))
c_eval('h1(?).FontSize = 14;',1:numel(h1))

%compact_panels(h2,0.0,00)
compact_panels(h2,0.005,00.005)
hb = findobj(gcf,'type','colorbar'); 
c_eval('hb(?).FontSize = 14;',1:numel(h2))
hb = hb(end:-1:1);
%
ihsub = [1:nt-1 nt+1:(2*nt-1) (2*nt-1)+1:(3*nt-1) (3*nt-1)+1:(4*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
%ihsub = [2:nt nt+2:(2*nt)];
ihsub = [2:nt nt+2:(2*nt) 2*nt+2:(3*nt) 3*nt+2:(4*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)


%% 2D distribution and stress contributions, MN, X times, with locations shown
ic = 3;
%t1 =  irf_time('2017-07-11T22:34:01.30Z','utc>EpochTT');
%t2 =  irf_time('2017-07-11T22:34:03.30Z','utc>EpochTT');
times = irf_time(['2017-07-11T22:34:01.30Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.03;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.4;
times = times(1):0.5:times(2);

times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z']);
%times = times([2 3 4])+0.25;
%times = times + 0.25;
times = times + 0.30;
%times = times + 0.06;
dt_dist = 2*0.062; % for two distributions

vint = [-Inf Inf];
vint = [-20000 20000];
vg = -60000:2000:60000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(3,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(1); 
h1.Position = position; 
%delete(h1_pos);
clear h2
for ip = 1:(nt*2)
  h2(ip) = subplot(3,nt,nt+ip);
end

hca = h1(1);
c_eval('mvaSe = mvaSe?;',ic)  
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaSe.xy*1e3,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaPe1.xy*1e3,mvaB1.abs},'comp')
hca.YLabel.String = 'S_{eLM} (pPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
h1.Position(4) = 0.16;
ppsum = [];
hp.EdgeColor = [0 0 0];
hp.FaceColor = [0.5 0.5 0.5];

times_exact = {};
isub = 1;
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca);
  hc_.Colorbar.YLabel.String = 'log_{10} f_e (s^2/m^5)';
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';

  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_distx = 0.030;
    B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(2)/b(1);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
    %B_
    irf_legend(hca,sprintf('B = %.1f nT', sqrt(sum(B_.^2))),[0.02 0.98],'color','k','fontsize',14)
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'stress');  
  hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';
  
  if 1
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    pptmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    ppsum(itime) = pptmp;
    
    if pptmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
  end

  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_dist = 0.030;
    B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(2)/b(1);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end
end

c_eval('hmark(?) = irf_pl_mark(h1(1),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]);',1:times.length)

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

h2(1).XLim = 0.99*vg([1 end])*1e-3;
h2(1).YLim = 0.99*vg([1 end])*1e-3;

c_eval('h2(?).XTick = -60:20:60; h2(?).YTick = -60:20:60;',1:numel(h2))

c_eval('h2(?).LineWidth = 1;',1:numel(h2))
c_eval('h2(?).FontSize = 14;',1:numel(h2))
c_eval('axis(h2(?),''square'');',1:numel(h2))

c_eval('h1(?).LineWidth = 1;',1:numel(h1))
c_eval('h1(?).FontSize = 14;',1:numel(h1))

%compact_panels(h2,0.0,00)
compact_panels(h2,0.005,00.005)
hb = findobj(gcf,'type','colorbar'); 
c_eval('hb(?).FontSize = 14;',1:numel(h2))
hb = hb(end:-1:1);
%
ihsub = [1:nt-1 nt+1:(2*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)

h1(1).Title.String = sprintf('MMS %g',ic);
h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
h1.Position(3) = h1.Position(3)*0.8;




if 1
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  tsPP = irf.ts_scalar(times,ppsum);
  hold(h1(1),'on')
  hpp = irf_plot(h1(1),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(1),'off')
end


%% PIC: example of plots, for explaining
%no02m = PIC('/Volumes/Fountain/cno062/data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists_new.h5');
figure(21);
twpe = 16000;
%zlim_ds = [-0.25 0.25];
%xlim_ds = [100.7 104.3];

%xpick = 102.4;
%xpick = 106.0;

zpick = 0;
for xpick = 102.4:0.2:106.0

  %ds = ds100.twpelim(twpe).xlim(xlim_ds).zlim(zlim_ds);
  ds = ds100.twpelim(twpe).xfind(xpick).zfind(zpick);
  clear h
  h(1) = subplot(2,1,1);
  h(2) = subplot(2,1,2);
  isub = 1;
  
  
  %ispecies = [2 4 6];
  ispecies = [4];
  sumdim = 3;
  fontsize = 14;
  
  hca = h(isub); isub = isub + 1;
  htmp = ds.plot_map(hca,ispecies,sumdim,'bline',no02m);
  hca.CLim = [0 0.002];
  
  
  hca = h(isub); isub = isub + 1;
  htmp = ds.plot_map(hca,ispecies,sumdim,'off-diag','bline',no02m);
  colormap(hca,pic_colors('blue_red'));
  hca.CLim = 0.04*[-1 1];
  
  ht = findobj(gcf,'type','text'); ht = ht(end:-1:1);
  delete(ht)
  
  if 1
    hi = findobj(hca.Children,'type','image');
    sumf = sum(hi.CData(:));
    if sumf > 0
      color = [1 0 0];
    else
      color = [0 0 1];
    end
    irf_legend(hca,sprintf('%.3f',sumf),[0.02 0.98],'color',color,'fontsize',fontsize)
    hca.Color = color;
  end
  
  %c_eval('h.ax(?).XGrid = ''off''; h.ax(?).YGrid = ''off'';',1:numel(h.ax))
  c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))
  
  
  for ip = 1:numel(h)
    axis(h(ip),'square');
    h(ip).XTick = -10:5:10;
    h(ip).YTick = -10:5:10;
    h(ip).XTickLabelRotation = 0;
  end
  c_eval('h(?).FontSize = fontsize;',1:numel(h))
  
  compact_panels(h,0.02)
  
  drawnow

  if 1 % print
    h(1).Title.String = sprintf('x = %.1f, z = %.1f',xpick,zpick);
    cn.print(sprintf('ex_fxy_x=%.1f_e4',xpick))
  
    
    h(1).XLabel.String = 'v_L';
    h(2).XLabel.String = 'v_L';
    h(1).YLabel.String = 'v_M';
    h(2).YLabel.String = 'v_M';
    cn.print(sprintf('ex_fLM_x=%.1f_e4',xpick))
  end
end

%% PIC: example of plots, for explaining, plot location of picked boxes
figure(19);
twpe = 16000;
%zlim_ds = [-0.25 0.25];
%xlim_ds = [100.7 104.3];

xpick = 102.4:0.2:106.0;
zpick = 0;

pic = no02m.twpelim(twpe).xlim(mean(no02m.xi)+[-4 4]).zlim([-2 2]);

h = pic.plot_map({'pexy'},'A',0.1,'sep');
h.CLim = 0.012*[-1 1];
colormap(pic_colors('blue_red'));
hold(h(1),'on')
[hd,hdl] = ds100.twpelim(twpe).xfind(xpick).zfind(zpick).plot_boxes(h(1));
hold(h(1),'off')
h(1).Position(2) = 0.18;
h(1).Position(4) = 0.7;
h(1).XLabel.String = 'L (d_i)';
h(1).YLabel.String = 'N (d_i)';
h(1).FontSize = 14;

hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
hb.YLabel.String = 'P^e_{LM}';

for ii = 1:numel(hdl)
  c_eval('hdl(?).LineWidth = 1;',1:numel(hdl))
  c_eval('hdl(?).LineWidth = 2;',ii)
  cn.print(sprintf('ex_locationboxes_x=%.1f',xpick(ii)))
end

%% PIC: example of plots, for explaining, combined with location boxes
%no02m = PIC('/Volumes/Fountain/cno062/data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists_new.h5');
figure(18);
twpe = 16000;
%zlim_ds = [-0.25 0.25];
%xlim_ds = [100.7 104.3];

%xpick = 102.4;
%xpick = 106.0;

doVertical = 1;
ispecies = [4 6];
zpick = 0;
sumdim = 3;
fontsize = 16;
pic = no02m.twpelim(twpe).xlim(mean(no02m.xi)+[-3 3]).zlim([-0.99 0.99]);

xpick_all = 102.4:0.2:106.0;
for xpick = xpick_all

  %ds = ds100.twpelim(twpe).xlim(xlim_ds).zlim(zlim_ds);
  ds = ds100.twpelim(twpe).xfind(xpick).zfind(zpick);
  clear h
  if not(doVertical)
    h(1) = subplot(1,4,[1 2]);
    h(2) = subplot(1,4,3);
    h(3) = subplot(1,4,4);
  
    h(1).Position(1) = 0.17;
    h(1).Position(3) = 0.3;
    
    %h(1).Position(1) = 0.08;
  else
    h(1) = subplot(3,4,[1 2 3 4]);
    h(2) = subplot(3,4,[5 6 9 10]);
    h(3) = subplot(3,4,[7 8 11 12]);
    h(3).Position(1) = h(3).Position(1) - 0.057;
    
    %h(1) = axes('Position',[0.1300    0.6593    0.7179    0.2157]);
    %h(2) = axes('Position',[0.1300    0.1100    0.3025    0.5154]);
    %h(3) = axes('Position',[0.5225    0.1100    0.3025    0.5154]);
    %h(1).Position(1) = 0.17;
    %h(1).Position(3) = 0.3;

  end
  isub = 1;

  % Location of boxes
  hca = h(isub); isub = isub + 1;
  hh = pic.plot_map(hca,{'pexy'},'A',0.1,'sep','smooth',3);
  hca.CLim = 0.012*[-1 1];
  colormap(hca,pic_colors('blue_red'));
  hold(hca,'on')
  %[hd,hdl] = ds100.twpelim(twpe).xfind(xpick).zfind(zpick).plot_boxes(hca);
  [hd,hdl] = ds100.twpelim(twpe).xfind(xpick).zfind(zpick).plot_boxes(hca);
  hold(hca,'off')
  %hca.Position(2) = 0.18;
  %hca.Position(4) = 0.7;
  hca.XLabel.String = 'L (d_i)';
  hca.YLabel.String = 'N (d_i)';
  hca.FontSize = 14;
  
  hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
  hb.YLabel.String = 'P^e_{LM}';
  
  
  
  hdl.LineWidth = 2;
  
  
  % f
  hca = h(isub); isub = isub + 1;
  htmp = ds.plot_map(hca,ispecies,sumdim,'bline',no02m);
  hca.CLim = [0 0.002];
  
  % Pressure integrand
  hca = h(isub); isub = isub + 1;
  htmp = ds.plot_map(hca,ispecies,sumdim,'off-diag',1/100,'bline',no02m);
  colormap(hca,pic_colors('blue_red'));
  hca.CLim = 0.04*[-1 1];
  
  ht = findobj(gcf,'type','text'); ht = ht(end:-1:1);
  delete(ht)
  
  if 1
    hi = findobj(hca.Children,'type','image');
    sumf = sum(hi.CData(:));
    if sumf > 0
      color = [1 0 0];
    else
      color = [0 0 1];
    end
    irf_legend(hca,sprintf('%.3f',sumf*1/100),[0.98 0.98],'color',color,'fontsize',fontsize)
    hca.Color = color;
  end
  
  %c_eval('h.ax(?).XGrid = ''off''; h.ax(?).YGrid = ''off'';',1:numel(h.ax))
  c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))
  
  
  for ip = 2:numel(h)
    axis(h(ip),'square');
    h(ip).XTick = -10:5:10;
    h(ip).YTick = -10:5:10;
    h(ip).XTickLabelRotation = 0;
  end
  c_eval('h(?).FontSize = fontsize;',1:numel(h))
  
  %compact_panels(h,0.02)
  
  drawnow


    h(2).XLabel.String = 'v_L';
    h(3).XLabel.String = 'v_L';
    h(2).YLabel.String = 'v_M';
    h(3).YLabel.String = 'v_M';

    if doVertical
      h(1).Position(2) = h(1).Position(2)-0.05;
      compact_panels(h(2:3),0,0.01)
      h(3).YLabel.String = [];
      h(3).YTickLabels = [];
    end
  if 1 % print
    %h(1).Title.String = sprintf('x = %.1f, z = %.1f',xpick,zpick);
  %  cn.print(sprintf('ex_fxy_x=%.1f_comb',xpick))
  
    
    cn.print(sprintf('ex_fLM_x=%.1f_comb_vertical',xpick))
  end
end

%% Map of off-diag PeLM
% for LM
xpicks = 102.4:0.4:105.0;
zpicks = -0.4:0.2:0;

% for MN
xpicks = 102.4:0.4:105.0;
zpicks = -0.2:0.1:0;

ds = ds100.twpelim(16000).xfind(xpicks).zfind(zpicks);


figure(103)
clear h
h(1) = subplot(1,1,1);
isub = 1;
hca = h(isub); isub = isub + 1;
hh = pic.plot_map(hca,{'peyz'},'A',0.1,'sep','smooth',3);
hca.CLim = 0.012*[-1 1];
colormap(hca,pic_colors('blue_red'));
hold(hca,'on')
%[hd,hdl] = ds100.twpelim(twpe).xfind(xpick).zfind(zpick).plot_boxes(hca);
[hd,hdl] = ds.plot_boxes(hca);
hold(hca,'off')
%hca.Position(2) = 0.18;
%hca.Position(4) = 0.7;
hca.XLabel.String = 'L (d_i)';
hca.YLabel.String = 'N (d_i)';
hca.FontSize = 14;

hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
hb.YLabel.String = 'P^e_{MN}';


%%
figure(104)
sumdim = 1;

h = findobj(gcf,'type','axes'); h = h(end:-1:1); delete(h);
h = ds.plot_map([2 4 6],sumdim,'bline',pic);
compact_panels(h.ax,0.01,0.01)
%c_eval('axis(h.ax(?),''equal'');',1:numel(h.ax))
c_eval('axis(h.ax(?),''square'');',1:numel(h.ax))
c_eval('h.ax(?).XGrid = ''off''; h.ax(?).YGrid = ''off'';',1:numel(h.ax))
hlinks = linkprop(h.ax,{'CLim'});

%hlinks.Targets(1).CLim = 0.04*[-1 1];
hlinks.Targets(1).CLim = 10e-3*[0 1];
c_eval('h.ax(?).XTick = -10:5:10;',1:numel(h.ax))
c_eval('h.ax(?).YTick = -10:5:10;',1:numel(h.ax))
c_eval('h.ax(?).XTickLabelRotation = 90;',1:numel(h.ax))

if sumdim == 3
  c_eval('h.ax(?).YLabel.String = ''v_M'';',1:numel(zpicks))
  c_eval('h.ax(?).XLabel.String = ''v_L'';',1:numel(zpicks):numel(h.ax))
  hb = findobj(gcf,'type','colorbar');
  hb.YLabel.String = 'f(v_M,v_L)';
elseif sumdim == 1
  c_eval('h.ax(?).YLabel.String = ''v_N'';',1:numel(zpicks))
  c_eval('h.ax(?).XLabel.String = ''v_M'';',1:numel(zpicks):numel(h.ax))
  hb = findobj(gcf,'type','colorbar');
  hb.YLabel.String = 'f(v_M,v_N)';


end
ht = findobj(gcf,'type','text'); ht = ht(end:-1:1);
delete(ht)

compact_panels(h.ax,0.003,0.003)

c_eval('h.ax(?).FontSize = 13;',1:numel(h.ax))


for ip = 1:numel(h.ax)
  hca = h.ax(ip);
  hi = findobj(hca.Children,'type','image');
  sumf = sum(hi.CData(:));
  if sumf > 0
    color = [1 0 0];
  else
    color = [0 0 1];
  end
  irf_legend(hca,sprintf('%.3f',sumf*1/100),[0.98 0.98],'color',color,'fontsize',fontsize)
  hca.Color = color;
  %colormap(hca,pic_colors('blue_red'))
  colormap(hca,pic_colors('candy4'))
  
end

%% Coordinate system
% illustrate_magnetic_reconnection

doVideo = 1;
doGif = 1;
fileName = 'illustration_magnetic_reconnection';

a = 5;
b = 1;
x = a*linspace(-10,10,200);
y = b*linspace(-10,10,100);
z = linspace(-10,10,5);
[X,Y] = meshgrid(x,y);
dx = x(2) - x(1);
dy = y(2) - y(1);
%dz = z(2) - z(1);
x_xline = x;
y_xline = x*b/a;

Ay = @(x,y) (x/a).^2 - (y/b).^2;
AY0 = Ay(X,Y);

%[FX,FY] = gradient(AY,dx,dy);
%Bx = -FX;
%By = FY;

colors = pic_colors('matlab');
colors = [colors; colors(end:-1:1,:)];

hca = subplot(1,1,1);
t = 0:30;
Astep = 20;
dA = Astep/numel(t);
AYlev0 = -100:Astep:(100 + Astep);

% Initiate

it = 1;
% Draw separatrix
if 0 % 2D
  plot(hca,x_xline,y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
  hold(hca,'on')
  plot(hca,x_xline,-y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
  % Draw field lines
  AY = AY0 - dA*t(it);
  S = contourcs(x,y,AY,AYlev0);
  for is = 1:numel(S)
    sx = interp1(1:numel(S(is).X),S(is).X,1:0.5:numel(S(is).X));
    sy = interp1(1:numel(S(is).Y),S(is).Y,1:0.5:numel(S(is).Y));%S(is).Y;
    plot(hca,sx,sy,'color',[0 0 0],'linewidth',2)
  end
elseif 1 % 3D
  plot3(hca,x_xline*0,x_xline,y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
  hold(hca,'on')
  plot3(hca,x_xline*0,x_xline,-y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
  % Draw field lines
  AY = AY0 - dA*t(it);
  S = contourcs(x,y,AY,AYlev0);
  for is = 1:numel(S)
    sx = interp1(1:numel(S(is).X),S(is).X,1:0.5:numel(S(is).X));
    sy = interp1(1:numel(S(is).Y),S(is).Y,1:0.5:numel(S(is).Y));%S(is).Y;
    plot3(hca,sx*0,sx,sy,'color',[0 0 0],'linewidth',2)
  end
end

pause(0.1)
drawnow
hold(hca,'off')
hca.Visible = 'off';
hca.Position = [0 0 1 1];

%% PIC: map of pxy
%no02m = PIC('/Volumes/Fountain/cno062/data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists_new.h5');
figure(21);
twpe = 16000;
%zlim_ds = [-0.25 0.25];
%xlim_ds = [100.7 104.3];
xpick = 100.8:0.4:104.2;
xpick = 100.0:0.4:105.0;
zpick = -0.3:0.1:0.3;

%ds = ds100.twpelim(twpe).xlim(xlim_ds).zlim(zlim_ds);
ds = ds100.twpelim(twpe).xfind(xpick).zfind(zpick);

%ispecies = [2 4 6];
ispecies = [6];
%h = ds.plot_map(ispecies,1,'off-diag','bline',no02m);
%h = ds.plot_map(ispecies,3,'bline',no02m);
h = ds.plot_map(ispecies,1,'ratio',[4 6],'bline',no02m);

compact_panels(h.ax,0.,0.)
%colormap(pic_colors('blue_red'));
c_eval('h.ax(?).XGrid = ''off''; h.ax(?).YGrid = ''off'';',1:numel(h.ax))




%
ht = findobj(gcf,'type','text'); ht = ht(end:-1:1);
delete(ht)
for ip = 1:numel(h.ax)
  axis(h.ax(ip),'square');
  h.ax(ip).XTick = -10:5:10;
  h.ax(ip).YTick = -10:5:10;
  h.ax(ip).XTickLabelRotation = 0;
  
  % sum of integrand
  if 0
    hi = findobj(h.ax(ip).Children,'type','image');
    sumf = sum(hi.CData(:));
    if sumf > 0
      color = [1 0 0];
    else
      color = [0 0 1];
    end
    irf_legend(h.ax(ip),sprintf('%.3f',sumf),[0.02 0.98],'color',color)
    h.ax(ip).Color = color;
  end
end

%% Plot location of boxes
figure(22);
pic = no02m.twpelim(twpe).xlim(mean(no02m.xi)+[-4 4]).zlim([-2 2]);

h = pic.plot_map({'pexy'});
h.CLim = 0.012*[-1 1];
colormap(pic_colors('blue_red'));
hold(h(1),'on')
ds.plot_boxes(h(1));
hold(h(1),'off')

%% PIC: peij
pic = no02m;

twpe = 16000;
xlim = mean(pic.xi) + 3*[-1 1];
zlim = 1.5*[-1 1];

varstrs = {'pexy','peyz'}';
cbarlabels = {'P^e_{LM}','P^e_{MN}'};
clims = {9.9e-3*[-1 1],9.9e-3*[-1 1]};
cmaps = {pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red')};

h = pic.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map(varstrs,'A',0.2,'sep','smooth',0,'clim',clims,'cbarlabels',cbarlabels,'cmap',cmaps);

c_eval('h(?).FontSize = 14;',1:numel(h))

c_eval('h(?).YLabel.String = ''N (d_i)'';',1:numel(h))
c_eval('h(?).XLabel.String = ''L (d_i)'';',1:numel(h))

%% PIC: peij, + gradients
pic = no02m;

twpe = 16000;
xlim = mean(pic.xi) + 3*[-1 1];
zlim = 1.5*[-1 1];

varstrs = {'pexy','peyz';'dxPexy','dzPezy'}';
cbarlabels = {'P^e_{LM}','P^e_{MN}';'\partial_LP^e_{LM}','\partial_NP^e_{MN}'}';
%cbarlabels = {'P^e_{LM}','P^e_{MN}';'\partial P^e_{LM}/\partial L','\partial P^e_{MN}/\partial N'}';

clims = {9.9e-3*[-1 1],9.9e-3*[-1 1];0.049*[-1 1],0.049*[-1 1]}';
cmaps = {pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red')};

h = pic.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map(varstrs,'A',0.2,'sep','smooth',0,'clim',clims,'cbarlabels',cbarlabels,'cmap',cmaps);

c_eval('axis(h(?),''equal'');',1:numel(h))
compact_panels(h,0.01,0.11)

hb = findobj(gcf,'type','colorbar'); ht = ht(end:-1:1);
c_eval('hb(?).FontSize = 16;',1:numel(hb))
c_eval('h(?).FontSize = 16;',1:numel(h))

c_eval('h(?).YLabel.String = ''N (d_i)'';',1:numel(h))
c_eval('h(?).XLabel.String = ''L (d_i)'';',1:numel(h))

% hl = findobj(gcf,'type','contour');
% c_eval('hl(?).LineWidth = 1;',1:numel(hl))

%% PIC: peij, + gradients divided by density
pic = no02m;

twpe = 16000;
xlim = mean(pic.xi) + 3*[-1 1];
%zlim = 1.5*[-1 1];
zlim = 0.99*[-1 1];

varstrs = {'pexy','peyz';'-dxPexy./ne','-dzPezy./ne'}';
cbarlabels = {'P^e_{LM}','P^e_{MN}';'-\partial_LP^e_{LM}/n','-\partial_NP^e_{MN}/n'}';
%cbarlabels = {'P^e_{LM}','P^e_{MN}';'\partial P^e_{LM}/\partial L','\partial P^e_{MN}/\partial N'}';

clims = {9.9e-3*[-1 1],9.9e-3*[-1 1];0.199*[-1 1],0.199*[-1 1]}';
cmaps = {pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red')};

h = pic.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map(varstrs,'A',0.2,'sep','smooth',3,'clim',clims,'cbarlabels',cbarlabels,'cmap',cmaps);

c_eval('axis(h(?),''equal'');',1:numel(h))
compact_panels(h,0.01,0.11)

hb = findobj(gcf,'type','colorbar'); ht = ht(end:-1:1);
c_eval('hb(?).FontSize = 16;',1:numel(hb))
c_eval('h(?).FontSize = 16;',1:numel(h))

c_eval('h(?).YLabel.String = ''N (d_i)'';',1:numel(h))
c_eval('h(?).XLabel.String = ''L (d_i)'';',1:numel(h))

%hl = findobj(gcf,'type','contour');
%c_eval('hl(?).LineWidth = 1;',1:numel(hl))
%c_eval('h(?).LineWidth = 1;',1:numel(h))
%c_eval('h(?).Position(2) = h(?).Position(2) + 0.02;',1:numel(h))

%% PIC: electron ohms law
pic = no02m;

twpe = 16000;
xlim = mean(pic.xi) + 3*[-1 1];
zlim = 0.99*[-1 1];

varstrs = {'Ey','vexBy','-(1/100)*vdvey','-dxPexy./ne','-dzPezy./ne'}';
cbarlabels = {'E','vxB','-m v\cdot\nabla v','-\partial_LP^e_{LM}/ne','-\partial_NP^e_{MN}/ne','sum'}';

varstrs = {'Ey+vexBy','-(1/100)*vdvey','-dxPexy./ne','-dzPezy./ne'}';
cbarlabels = {'E+vxB','-m v\cdot\nabla v','-\partial_LP^e_{LM}/ne','-\partial_NP^e_{MN}/ne','sum'}';

%varstrs = {'Ey+vexBy','-(1/100)*vdvey-dxPexy./ne-dzPezy./ne'}';
%cbarlabels = {'E+vxB','-v\cdot\nabla v - \partial_LP^e_{LM}/n - \partial_NP^e_{MN}/n','sum'}';


clims = {0.199*[-1 1],0.199*[-1 1],0.199*[-1 1],0.199*[-1 1],0.199*[-1 1],0.199*[-1 1]}';
cmapbr = pic_colors('blue_red');
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr};

h = pic.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map(varstrs,'A',0.2,'sep','smooth',3,'clim',clims,'cbarlabels',cbarlabels,'cmap',cmaps);

%c_eval('axis(h(?),''equal'');',1:numel(h))
%compact_panels(h,0.01,0.11)

hb = findobj(gcf,'type','colorbar'); ht = ht(end:-1:1);
c_eval('hb(?).FontSize = 16;',1:numel(hb))
c_eval('h(?).FontSize = 16;',1:numel(h))

c_eval('h(?).YLabel.String = ''N (d_i)'';',1:numel(h))
c_eval('h(?).XLabel.String = ''L (d_i)'';',1:numel(h))

c_eval('h(?).Position(3) = 0.7;',1:numel(h))

if 0

  %%
  hl = findobj(gcf,'type','contour');
  c_eval('hl(?).LineWidth = 0.5;',1:numel(hl))
  c_eval('h(?).LineWidth = 0.5;',1:numel(h))
end

%% PIC: electron anisotropy
%no02m = PIC('/Volumes/Fountain/cno062/data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists_new.h5');
figure(22);
twpe = 16000;

pic = no02m.twpelim(twpe).xlim(mean(no02m.xi)+[-4 4]).zlim([-2 2]);

varstrs = {'log10(tepar./teperp)','tepar./teperp-1'}';
varstrs = {'log10(tepar./teperp)'}';
cbarlabels = {'log_{10}T^e_{||}/T^e_{\perp}'}';
h = pic.plot_map(varstrs,'A',0.2,'sep','cbarlabels',cbarlabels);
colormap(pic_colors('blue_red'));

c_eval('h(?).FontSize = 16;',1:numel(h))

c_eval('h(?).YLabel.String = ''N (d_i)'';',1:numel(h))
c_eval('h(?).XLabel.String = ''L (d_i)'';',1:numel(h))

h(1).CLim = 1.5*[-1 1];

if 1
  %%
  
  hca = h(1);
  xpicks = 102.4;
  zpicks = -0.4;
  ds = ds100.twpelim(16000).xfind(xpicks).zfind(zpicks);
  hold(hca,'on')
  %[hd,hdl] = ds100.twpelim(twpe).xfind(xpick).zfind(zpick).plot_boxes(hca);
  [hd,hdl] = ds.plot_boxes(hca);
  hold(hca,'off')
end

%% Compare pressure to stress tensor

h = irf_plot(8);

hca = irf_panel('ve');
c_eval('irf_plot(hca,{mvaVe?.x,mvaVe?.y,mvaVe?.z},''comp'');',ic)
hca.YLabel.String = 'v_e (km/s)';

hca = irf_panel('Te anis');
c_eval('irf_plot(hca,{Temprat?},''comp'');',ic)
hca.YLabel.String = 'T^e_{||}/T^e_{\perp}';

for comp = ["xx","yy","zz","xy","xz","yz"]
  hca = irf_panel(char(comp));
  c_eval('irf_plot(hca,{mvaPe?.(comp),mvaSe?.(comp)},''comp'')',ic)
  hca.YLabel.String = comp;
  irf_legend(hca,{'P','S'}',[1.02,0.95])
end
