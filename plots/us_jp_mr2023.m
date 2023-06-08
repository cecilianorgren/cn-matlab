%% US-Japan
% Load data
ic = 1;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854

% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
mms.db_init('local_file_db','/Users/cno062/Data/MMS');
db_info = datastore('mms_db');

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
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',ic);

% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);

% Particle moments
% Skymap distributions
disp('Loading skymaps...')
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic)
c_eval('iPDist?_onecount = iPDist?; iPDist?_onecount.data = (iPDist?_onecount.data./iPDistErr?.data).^2;',ic)

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
% Rotate things into new coordinate system 
L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9239];
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

c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;',ic)


c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';',ic)

%% Figure: Overview 1
ic = 1;

tint_edr = irf.tint('2017-07-11T22:33:58.00Z/2017-07-11T22:34:10.00Z'); %20151112071854

npanels = 4;
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
  hca = irf_panel('Pe off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPe?.xy.tlim(tint)*1e3,mvaPe?.xz.tlim(tint)*1e3,mvaPe?.yz.tlim(tint)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'P_e','(pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
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
ic = 1;
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
times = times + 0.25;

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
irf_plot(hca,mvaPe1.xy*1e3,'color','k','linewidth',1)
hca.YLabel.String = 'P_{eLM} (pPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
h1.Position(4) = 0.16;

times_exact = {};
isub = 1;
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
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
    xlim = hca.XLim;
    tlim = hca.YLim;
    hold(hca,'on')
    dt_dist = 0.030;
    B = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B = mean(B.data,1);
    b = B/norm(B);
    k = b(2)/b(1);
    if k > 1
      plot(hca,xlim,xLim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',mvaVe1.x.resample(dist),mvaVe1.y.resample(dist));  
  hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';

  if 1 % plot B direction
    xlim = hca.XLim;
    tlim = hca.YLim;
    hold(hca,'on')
    dt_dist = 0.030;
    B = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B = mean(B.data,1);
    b = B/norm(B);
    k = b(2)/b(1);
    if k > 1
      plot(hca,xlim,xLim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end
end

c_eval('hmark(?) = irf_pl_mark(h1(1),times_exact{?}.epochUnix + 0.5*0.03*[-1 1],mms_colors(''1''));',1:times.length)

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