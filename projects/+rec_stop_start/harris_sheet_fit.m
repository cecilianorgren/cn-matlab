%% Specify time and spacecraft
units = irf_units;
irf.log('critical')
ic = 1:4;

mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
tint_all = irf.tint('2017-01-01T00:00:00.00Z/2018-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

fileId = '20170725220853';
iFile = find(cellfun(@(s) contains(s,fileId),{files.name}));
tint = [files(iFile-1).start files(iFile).stop] + [1 -1];

tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');

% Event path
eventPath = ['/Users/' localuser '/Research/Events/mms_' fileId '/']; % for saving figures
mkdir(eventPath)

%% Load data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('gsmE? = c_coord_trans(''GSE'',''GSM'',gseE?);',ic)
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('gsmVe? = c_coord_trans(''GSE'',''GSM'',gseVe?);',ic)
c_eval('gsmVi? = c_coord_trans(''GSE'',''GSM'',gseVi?);',ic)
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4)
c_eval('gsmVExB? = cross(gsmE?.resample(gsmB?.time),gsmB?)/gsmB?.abs/gsmB?.abs*1e3; gsmVExB?.units = '''';',ic) % km/s
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('[gseVi?par,gseVi?perp] = irf_dec_parperp(gseB?,gseVi?); gseVi?par.name = ''Vi par''; gseVi?perp.name = ''Vi perp'';',ic)

gseBav = (gseB1 + gseB4.resample(gseB1) + gseB3.resample(gseB1) + gseB4.resample(gseB1))/4;

%% Current, curlometer from 4 sc
if all(ic == [1:4])
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurl.data = gseJcurl.data; Jcurl.units = 'nAm^{-2}';
gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';
end
% Currents from moments, use ne also for Ji 
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);

gseJav = (gseJ1 + gseJ2.resample(gseJ1) + gseJ3.resample(gseJ1) + gseJ4.resample(gseJ1))/4;
gseJiav = (gseJi1 + gseJi2.resample(gseJi1) + gseJi3.resample(gseJi1) + gseJi4.resample(gseJi1))/4;
gseJeav = (gseJe1 + gseJe2.resample(gseJe1) + gseJe3.resample(gseJe1) + gseJe4.resample(gseJe1))/4;

%% Figure, 4 sc
tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');
%tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:07:30.00Z');
%tint_harris = irf.tint('2017-07-25T22:07:30.00Z/2017-07-25T22:09:00.00Z');
%tint_harris = irf.tint('2017-07-25T22:09:00.00Z/2017-07-25T22:09:30.00Z');
units = irf_units;
% Parameters for Harris fit. Want an expression Bx(Jy)
syms z z0 L B0
Bx = B0*tanh((z-z0)/L);
Jy = gradient(Bx,z)/units.mu0; % -(B0*(tanh((z - z0)/L)^2 - 1))/L
JxBz = -Jy*Bx;

% Jy = -(B0*(tanh((z - z0)/L)^2 - 1))/L 
%    = -(B0*((Bx/B0)^2 - 1))/L
%     
% L  = -(B0*((Bx/B0)^2 - 1))/Jy
mf_L = @(Bx,B0,Jy) -(B0*((Bx/B0)^2 - 1))/Jy;

mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);
mf_JxBz = matlabFunction(JxBz);

mf_L_B = solve(Bx,L);
mf_L_J = solve(Jy,L);

if 0 % Illustration of Harris sheet
  %%
  h = setup_subplots(2,3);
  z0_ = 00e3;
  ii = 0;
  legs = {};
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500]*1e3  
      ii = ii + 1;      
      isub = 1;
      legs{ii} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,50);
      
      hca = h(isub); isub = isub + 1;
      plot(hca,z_*1e-3,mf_Bx(B0_,L_,z_,z0_)*1e9)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'B_x (nT)';
      
      hca = h(isub); isub = isub + 1;
      plot(hca,z_*1e-3,mf_Jy(B0_,L_,z_,z0_)*1e9)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'J_y (nA/m^2)';
      
      hca = h(isub); isub = isub + 1;
      plot(hca,z_*1e-3,mf_JxBz(B0_,L_,z_,z0_)*1e9)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'JxB_z (T nA/m^2)';
      
      hca = h(isub); isub = isub + 1;
      plot(hca,mf_Bx(B0_,L_,z_,z0_)*1e9,mf_Jy(B0_,L_,z_,z0_)*1e9)
      hca.XLabel.String = 'B_x (nT)';
      hca.YLabel.String = 'J_y (nA/m^2)';
      
      hca = h(isub); isub = isub + 1;
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9)
      hca.XLabel.String = 'B_x^2 (nT)';
      hca.YLabel.String = 'J_y (nA/m^2)';
         
      hca = h(isub); isub = isub + 1;
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxBz(B0_,L_,z_,z0_)*1e9)
      hca.XLabel.String = 'B_x^2 (nT)';
      hca.YLabel.String = 'JxB_z (T nA/m^2)';
      
      if ii == 1, c_eval('hold(h(?),''on'');',1:numel(h)); end
    end
  end  
  c_eval('hold(h(?),''off'');',1:4)
  c_eval('h(?).XGrid = ''on'';',1:4)
  c_eval('h(?).YGrid = ''on'';',1:4)
  %h(3).XLim(1) = 0;
  legend(h(1),legs,'box','off','location','southeast')
end
ic = 2;

npanels = 4;
nrows = 2;
ncols = 1;
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');
iisub = 0;
cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 1 % J, Jeav, Jiav, Jav ,curl
  for comp = ['x','y','z']      
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k')
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',12);  
  end
end
if 0 % J, mms1-4,curl
  for comp = ['x','y','z']  
    species = 'i';    
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    %irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    eval(sprintf('irf_plot(hca,{gseJ%s1.(comp),gseJ%s2.(comp),gseJ%s3.(comp),gseJ%s4.(comp),gseJcurl.(comp)},''comp'');',species,species,species,species))
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s,%s}^{mms1}',species,comp),sprintf('J_{%s,%s}^{mms2}',species,comp),sprintf('J_{%s,%s}^{mms3}',species,comp),sprintf('J_{%s,%s}^{mms4}',species,comp),'curl'},[0.02 0.99],'fontsize',12);  
  end
end

irf_zoom(h1,'x',tint)
irf_zoom(h1(zoomy),'y')
irf_plot_axis_align
irf_pl_mark(h1(1),tint_harris)

% Non-TSeries panels.
isub = 1;
if 0 % Bx,Jy
  hca = h2(isub); isub = isub + 1;
  plot(hca,abs(gseBav.tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'.')  
end

for species = {'','i','e'}
  if 0 % Jy,Jcurl
    hca = h2(isub); isub = isub + 1;
    eval(sprintf('xx = gseJ%sav.tlim(tint).y.data;',species{1}))
    yy = gseJcurl.resample(gseJav).tlim(tint).y.data;
    irem = find(xx == 0);
    xx(irem) = [];
    yy(irem) = [];
    p = polyfit(xx,yy,1);
    J_edges1 = linspace(-10,10,41);
    J_edges2 = linspace(-10,10,39);
    [N edges mid loc] = histcn([xx(:) yy(:)],J_edges1,J_edges2);
    N(N==0) = NaN;
    pcolor(hca,mid{:},log10(N)')
    shading(hca,'flat')
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'log_{10}(counts)';
    colormap(pic_colors('candy4'))
    %plot(hca,xx,yy,'.')
    %histcn_plot(hca,xx,yy)
    hca.XLabel.String = sprintf('<J_{%sy}^{FPI}>',species{1});
    hca.YLabel.String = 'J_y^{curl}';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hold(hca,'on')
    xlim = hca.XLim;
    ylim = hca.YLim;
    plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
    plot(hca,[-100 100],p(2)+[-100 100]*p(1))
    hca.XLim = xlim;
    hca.YLim = ylim;
    hold(hca,'off')
    %hca.XLim = 40*[-1 1];
    %hca.YLim = 40*[-1 1]; 
    hca.XGrid = 'on';
    hca.YGrid = 'on'; 
  end
end
if 0 % Jey,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJiav.tlim(tint).y.data;
  yy = gseJcurl.resample(gseJiav).tlim(tint).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{i,y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 0 % Jiy,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJeav.tlim(tint).y.data;
  yy = gseJcurl.resample(gseJeav).tlim(tint).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{e,y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 1 % Bx,Jiy
  hca = h2(isub); isub = isub + 1;
  %plot(hca,abs(gseBav.resample(gseJi1).tlim(tint_harris).x.data),gseJi1.tlim(tint_harris).y.data,'.')
  %plot(hca,abs(gseBav.resample(gseJcurl).tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'k.')
  %ts_yy = gseJiav; y_str = '<J_{iy}^{FPI}>';
  %ts_yy = gseJcurl.resample(timeline); y_str = '<J_{y}^{curl}>';
  ts_yy = gseJcurl; y_str = '<J_{y}^{curl}>';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  %irem = find(xx == 0);
  %xx(irem) = [];
  %yy(irem) = [];
  %p = polyfit(xx,yy,1);
  J_edges = linspace(-2,12,70);
  B_edges = linspace(0,24,71);
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(counts)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   %% 
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end
if 1 % Bx,Jiy, time
  hca = h2(isub); isub = isub + 1;
  %plot(hca,abs(gseBav.resample(gseJi1).tlim(tint_harris).x.data),gseJi1.tlim(tint_harris).y.data,'.')
  %plot(hca,abs(gseBav.resample(gseJcurl).tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'k.')
  %ts_yy = gseJiav; y_str = '<J_{iy}^{FPI}>';
  %ts_yy = gseJcurl.resample(timeline); y_str = '<J_{y}^{curl}>';
  ts_yy = gseJcurl; y_str = '<J_{y}^{curl}>';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  
  scatter(hca,xx.^2,yy,1,tt)  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('time since %s',ts_yy.resample(ts_yy).tlim(tint_harris).time(1).utc);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   %% 
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end
if 0 % Bx,JxBz, binned
  hca = h2(isub); isub = isub + 1;
  %plot(hca,abs(gseBav.resample(gseJi1).tlim(tint_harris).x.data),gseJi1.tlim(tint_harris).y.data,'.')
  %plot(hca,abs(gseBav.resample(gseJcurl).tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'k.')
  %ts_yy = gseJiav; y_str = '<J_{iy}^{FPI}>';
  %ts_yy = gseJcurl.resample(timeline); y_str = '<J_{y}^{curl}>';
  ts_yy = gseJxB; y_str = '-J_{y}xB_x';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  %irem = find(xx == 0);
  %xx(irem) = [];
  %yy(irem) = [];
  %p = polyfit(xx,yy,1);
  JxB_edges = linspace(-2,12,70);
  B_edges = linspace(0,24,71);
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},log10(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(counts)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
    
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end
if 0 % Bx,JxBz, scatter
  hca = h2(isub); isub = isub + 1;
  
  ts_yy = gseJxB; y_str = '-J_{y}xB_x';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data*1e18;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  cc = ts_yy.tlim(tint_harris).time - ts_yy.tlim(tint_harris).time(1);
  scatter(hca,xx.^2,yy,1,cc)
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'time since tt';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
    
  if 0
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  end
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end

%% Figure, 4 sc
tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');
%tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:07:30.00Z');
%tint_harris = irf.tint('2017-07-25T22:07:30.00Z/2017-07-25T22:09:00.00Z');
%tint_harris = irf.tint('2017-07-25T22:09:00.00Z/2017-07-25T22:09:30.00Z');
untis = irf_units;
% Parameters for Harris fit. Want an expression Bx(Jy)
syms z z0 L B0
Bx = B0*tanh((z-z0)/L);
Jy = gradient(Bx,z)/units.mu0; % -(B0*(tanh((z - z0)/L)^2 - 1))/L

mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);

if 0 % Illustration of Harris sheet
  %%
  h = setup_subplots(2,2);
  z0_ = 00e3;
  ii = 0;
  legs = {};
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500]*1e3  
      ii = ii + 1;      
      legs{ii} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);
      hca = h(1);
      plot(hca,z_*1e-3,mf_Bx(B0_,L_,z_,z0_)*1e9)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'B_x (nT)';
      hca = h(2);
      plot(hca,z_*1e-3,mf_Jy(B0_,L_,z_,z0_)*1e9)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'J_y (nA/m^2)';
      hca = h(3);
      plot(hca,mf_Bx(B0_,L_,z_,z0_)*1e9,mf_Jy(B0_,L_,z_,z0_)*1e9)
      hca.XLabel.String = 'B_x (nT)';
      hca.YLabel.String = 'J_y (nA/m^2)';
      hca = h(4);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9)
      hca.XLabel.String = 'B_x^2 (nT)';
      hca.YLabel.String = 'J_y (nA/m^2)';
      
      if ii == 1, c_eval('hold(h(?),''on'');',1:4); end
    end
  end  
  c_eval('hold(h(?),''off'');',1:4)
  c_eval('h(?).XGrid = ''on'';',1:4)
  c_eval('h(?).YGrid = ''on'';',1:4)
  %h(3).XLim(1) = 0;
  legend(h(1),legs,'box','off','location','southeast')
end
ic = 2;

npanels = 4;
nrows = 2;
ncols = 1;
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.7,'vertical');
iisub = 0;
cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 0 % J, Jeav, Jiav, Jav ,curl
  for comp = ['x','y','z']
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k')
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',12);  
  end
end
if 1 % Jcurl, Jfpiav , Jcurl-Jfpiav
  for comp = ['x','y','z']
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('14b'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)-gseJav.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k')
    else
      irf_plot(hca,{gseJav.(comp),gseJcurl.(comp).resample(gseJav),gseJcurl.(comp).resample(gseJav)-gseJav.(comp)},'comp');
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('14b'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp),sprintf('J_{%s}^{curl}-J_{%s}^{FPI}',comp,comp)},[0.02 0.99],'fontsize',12);  
  end
end
if 0 % J, mms1-4,curl
  for comp = ['x','y','z']  
    species = 'i';    
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    %irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    eval(sprintf('irf_plot(hca,{gseJ%s1.(comp),gseJ%s2.(comp),gseJ%s3.(comp),gseJ%s4.(comp),gseJcurl.(comp)},''comp'');',species,species,species,species))
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s,%s}^{mms1}',species,comp),sprintf('J_{%s,%s}^{mms2}',species,comp),sprintf('J_{%s,%s}^{mms3}',species,comp),sprintf('J_{%s,%s}^{mms4}',species,comp),'curl'},[0.02 0.99],'fontsize',12);  
  end
end

irf_zoom(h1,'x',tint)
irf_zoom(h1(zoomy),'y')
irf_plot_axis_align
irf_pl_mark(h1(1),tint_harris)

% Non-TSeries panels.
isub = 1;
if 0 % Bx,Jy
  hca = h2(isub); isub = isub + 1;
  plot(hca,abs(gseBav.tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'.')  
end

for species = {''}
  if 1 % Jy,Jcurl
    hca = h2(isub); isub = isub + 1;
    eval(sprintf('xx = gseJ%sav.tlim(tint).y.data;',species{1}))
    yy = gseJcurl.resample(gseJav).tlim(tint).y.data;
    irem = find(xx == 0);
    xx(irem) = [];
    yy(irem) = [];
    p = polyfit(xx,yy,1);
    J_edges1 = linspace(-10,10,41);
    J_edges2 = linspace(-10,10,39);
    [N edges mid loc] = histcn([xx(:) yy(:)],J_edges1,J_edges2);
    N(N==0) = NaN;
    pcolor(hca,mid{:},log10(N)')
    shading(hca,'flat')
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'log_{10}(counts)';
    colormap(pic_colors('candy4'))
    %plot(hca,xx,yy,'.')
    %histcn_plot(hca,xx,yy)
    hca.XLabel.String = sprintf('<J_{%sy}^{FPI}>',species{1});
    hca.YLabel.String = 'J_y^{curl}';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hold(hca,'on')
    xlim = hca.XLim;
    ylim = hca.YLim;
    plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
    plot(hca,[-100 100],p(2)+[-100 100]*p(1))
    hca.XLim = xlim;
    hca.YLim = ylim;
    hold(hca,'off')
    %hca.XLim = 40*[-1 1];
    %hca.YLim = 40*[-1 1]; 
    hca.XGrid = 'on';
    hca.YGrid = 'on'; 
    axis(hca,'square')
  end
end
if 0 % Jey,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJiav.tlim(tint).y.data;
  yy = gseJcurl.resample(gseJiav).tlim(tint).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{i,y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 0 % Jiy,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJeav.tlim(tint).y.data;
  yy = gseJcurl.resample(gseJeav).tlim(tint).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{e,y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 1 % Bx,Jiy
  hca = h2(isub); isub = isub + 1;
  %plot(hca,abs(gseBav.resample(gseJi1).tlim(tint_harris).x.data),gseJi1.tlim(tint_harris).y.data,'.')
  %plot(hca,abs(gseBav.resample(gseJcurl).tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'k.')
  ts_yy = gseJiav; y_str = '<J_{iy}^{FPI}>';
  %ts_yy = gseJcurl.resample(timeline); y_str = '<J_{y}^{curl}>';
  %ts_yy = gseJcurl; y_str = '<J_{y}^{curl}>';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  %irem = find(xx == 0);
  %xx(irem) = [];
  %yy(irem) = [];
  %p = polyfit(xx,yy,1);
  J_edges = linspace(-2,12,70);
  B_edges = linspace(0,24,71);
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},log10(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(counts)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
    
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end

c_eval('h1(?).YLim = [-14 14];',2:4)

%% Figure, compare v to vExB
dt_resample = 0.5;
%timeline = tint(1):dt_resample:tint(2);
timeline = gseVi2.time; dt_resample = timeline(2) - timeline(1);
doResample = 1;


npanels = 3;
nrows = 3;
ncols = 2;
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');

ic = 2;
if 1 % Vi, Ve, VExB
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(sprintf('v%s',comp));
    set(hca,'ColorOrder',mms_colors('143'))
    
    c_eval('irf_plot(hca,{gseVe?perp.(comp).resample(timeline),gseVi?perp.(comp).resample(timeline),gseVExB?.(comp).resample(timeline)},''comp'');',ic)
    irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k')
    
    hca.YLabel.String = {['v_{\perp' comp '}'],'(km/s)'};
    set(hca,'ColorOrder',mms_colors('143'))
    irf_legend(hca,{'v_e','v_i','v_{ExB}'}',[1.02 0.99],'fontsize',12);   
    hca.YLim = [-1200 1200];
  end
end

irf_zoom(h1,'x',tint)
%irf_zoom(h1(zoomy),'y')
irf_plot_axis_align
c_eval('irf_pl_mark(h1(?),tint_harris)',1:numel(h))

isub = 1;
for species = ['i','e']
  for comp = ['x','y','z']
    hca = h2(isub); isub = isub + 1;
    c_eval('xx = gseVExB?.tlim(tint_harris).(comp).resample(timeline).data;',ic)
    c_eval(sprintf('yy = gseV%s?.tlim(tint_harris).(comp).resample(timeline).data;',species),ic)    
    irem = find(xx == 0);
    xx(irem) = [];
    yy(irem) = [];
    p = polyfit(xx,yy,1);
    J_edges1 = linspace(-900,900,41);
    J_edges2 = linspace(-900,900,39);
    [N edges mid loc] = histcn([xx(:) yy(:)],J_edges1,J_edges2);
    N(N==0) = NaN;
    pcolor(hca,mid{:},log10(N)')
    shading(hca,'flat')
    %hcb = colorbar('peer',hca);
    %hcb.YLabel.String = 'log_{10}(counts)';
    colormap(pic_colors('candy4'))
    %plot(hca,xx,yy,'.')
    %histcn_plot(hca,xx,yy)
    %hca.XLabel.String = ['v_{ExB ' comp '}'];
    hca.XLabel.String = ['v_{ExB}'];
    hca.YLabel.String = ['v_{' species ' ' comp '}'];
    irf_legend(hca,{['v_{' species ' ' comp '}']},[0.02 0.98]);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hold(hca,'on')
    xlim = hca.XLim;
    ylim = hca.YLim;
    plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
    plot(hca,[-100 100],p(2)+[-100 100]*p(1))
    hca.XLim = xlim;
    hca.YLim = ylim;
    hold(hca,'off')
    %hca.XLim = 40*[-1 1];
    %hca.YLim = 40*[-1 1]; 
    hca.XGrid = 'on';
    hca.YGrid = 'on'; 
    axis(hca,'square')  
  end
end
compact_panels(h2)

%% hodogram plot of Bav v Jcurl
timeline = tint_harris(1):0.01:tint_harris(2);
timeline = gseBav.tlim(tint_harris).time;

J_edges = linspace(-2,12,70);
B_edges = linspace(0,24,71);
%[BB,JJ] = 
B0 = 21; % nT
L = mf_L(B_edges,B0,J_edges);


h = setup_subplots(3,1); isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,gseBav.x.tlim(tint_harris).data.^2,gseJcurl.resample(gseBav).y.tlim(tint_harris).data);
hca.XLabel.String = 'B_x^2 [(nT)^2]';
hca.XLabel.String = 'B_x^2 [(nT)^2]';
hca.YLabel.String = 'J_{y}^2 [nA/m^2]';

hca = h(isub); isub = isub + 1;
hs = scatter(hca,gseBav.x.resample(timeline).data.^2,gseJcurl.resample(gseBav).y.resample(timeline).data,1,timeline-timeline(1));
hca.XLabel.String = 'B_x^2 [(nT)^2]';
hca.YLabel.String = 'J_{y}^2 [nA/m^2]';
hs.Marker = 'o';
hs.MarkerFaceAlpha = 1;
hs.MarkerEdgeAlpha = 1;
hs.MarkerFaceColor = hs.MarkerEdgeColor;
hcb = colorbar('peer',hca);
hcb.YLabel.String = sprintf('time since %s',timeline(1).utc);
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.YLim(2) = 20;

hca = h(isub); isub = isub + 1;
scatter(hca,gseBav.abs.tlim(tint_harris).data.^2,gseJcurlperp.resample(gseBav).y.tlim(tint_harris).data,1,timeline-timeline(1))
hca.XLabel.String = 'B^2 [(nT)^2]';
hca.YLabel.String = 'J_{\perp}^2 [nA/m^2]';
hs.Marker = 'o';
hs.MarkerFaceAlpha = 1;
hs.MarkerEdgeAlpha = 1;
hs.MarkerFaceColor = hs.MarkerEdgeColor;
hcb = colorbar('peer',hca);
hcb.YLabel.String = sprintf('time since %s',timeline(1).utc);
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.YLim(2) = 20;



%% From simulation
units = irf_units;
if 0
  %%
  no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
  Bx_orig = no02m.twpelim(24000).Bx;
  Jy_orig = no02m.twpelim(24000).Jy;
  jiy_orig = no02m.twpelim(24000).jiy;
  jey_orig = no02m.twpelim(24000).jey;
end
xlim = [72 93]; %xlim = [102.5 110]; %xlim = [no02m.xi(1) no02m.xi(end)];
xlim = [72 133];
xlim = [110 161];

zlim = [-10 10];
Bx = Bx_orig;
Bx(or(no02m.xi<xlim(1),no02m.xi>xlim(2)),:) = [];
Bx(:,or(no02m.zi<zlim(1),no02m.zi>zlim(2))) = [];

Jy = Jy_orig; Jy(abs(Jy)<0.001) = NaN;
Jy(or(no02m.xi<xlim(1),cno02m.xi>xlim(2)),:) = [];
Jy(:,or(no02m.zi<zlim(1),no02m.zi>zlim(2))) = [];

jiy = jiy_orig; 
jiy(or(no02m.xi<xlim(1),no02m.xi>xlim(2)),:) = [];
jiy(:,or(no02m.zi<zlim(1),no02m.zi>zlim(2))) = [];

jey = jey_orig; 
jey(or(no02m.xi<xlim(1),no02m.xi>xlim(2)),:) = [];
jey(:,or(no02m.zi<zlim(1),no02m.zi>zlim(2))) = [];

xi = no02m.xi;
zi = no02m.zi;
xi(or(no02m.xi<xlim(1),no02m.xi>xlim(2))) = [];
zi(or(no02m.zi<zlim(1),no02m.zi>zlim(2))) = [];


syms sz sz0 sL sB0
sBx = sB0*tanh((sz-sz0)/sL);
sJy = gradient(sBx,sz); % -(B0*(tanh((z - z0)/L)^2 - 1))/L

mf_Bx = matlabFunction(sBx);
mf_Jy = matlabFunction(sJy);


jielim = 0.3*[-1 1];
jlim = [-0.2 1];
blim = max(Bx_orig(:))*[-1 1];

h = setup_subplots(2,2,'vertical');
isub = 1;

if 1 % Bx
  hca = h(isub); isub = isub + 1;
  pcolor(hca,xi,zi,Bx')
  shading(hca,'flat')
  hcb = colorbar('peer',hca); 
  hcb.YLabel.String = 'B_x';
  hcb.YLim = max(blim)*[-1 1];
  hca.CLim = blim;
  colormap(hca,pic_colors('blue_red'))
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'z/d_i';
end
if 1 % Jy
  hca = h(isub); isub = isub + 1;
  pcolor(hca,xi,zi,Jy')
  shading(hca,'flat')
  hcb = colorbar('peer',hca); 
  hcb.YLabel.String = 'J_y';
  hcb.YLim = jlim;
  hca.CLim = max(jlim)*[-1 1];
  colormap(hca,pic_colors('blue_gray_red'))
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'z/d_i';
end
if 1 % (Bx^2,Jy)
  hca = h(isub); isub = isub + 1;  
  [N edges mid loc] = histcn([Bx(:).^2 Jy(:)],linspace(0,blim(2).^2,100),linspace(jlim(1),jlim(2),100));
  N(N==0) = NaN;
  pcolor(hca,mid{1:2},log10(N(:,:,ceil(end/2)))')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(counts)';
  colormap(hca,pic_colors('candy4'))
  hca.XLabel.String = 'B_x^2';
  hca.YLabel.String = 'J_y';
  hca.XGrid = 'on';
  hca.YGrid = 'on';  

  legs = {};
  hl = gobjects(0);
  hold(hca,'on');
  B0 = 0.28; L = 1.6; z0 = 0; zvec = linspace(0,3*L,100); legs{end+1} = sprintf('B_0 = %g , L = %g',B0,L);  
  hl(end+1) = plot(hca,mf_Bx(B0,L,zvec,z0).^2,mf_Jy(B0,L,zvec,z0));
  plot(hca,mf_Bx(B0,L,zvec,z0).^2,mf_Jy(B0,L,zvec,z0),'k:');
  
 
  B0 = 0.35; L = 0.5; z0 = 0; zvec = linspace(0,3*L,100); legs{end+1} = sprintf('B_0 = %g , L = %g',B0,L);  
  hl(end+1) = plot(hca,mf_Bx(B0,L,zvec,z0).^2,mf_Jy(B0,L,zvec,z0));
  plot(hca,mf_Bx(B0,L,zvec,z0).^2,mf_Jy(B0,L,zvec,z0),'k:');
  
  B0 = 1.05; L = 2; z0 = 0; zvec = linspace(0,3*L,100); legs{end+1} = sprintf('B_0 = %g , L = %g',B0,L);  
  hl(end+1) = plot(hca,mf_Bx(B0,L,zvec,z0).^2,mf_Jy(B0,L,zvec,z0));
  plot(hca,mf_Bx(B0,L,zvec,z0).^2,mf_Jy(B0,L,zvec,z0),'k:');
  
  B0 = 0.50; L = 3; z0 = 0; zvec = linspace(0,3*L,100); legs{end+1} = sprintf('B_0 = %g , L = %g',B0,L);
  hl(end+1) = plot(hca,mf_Bx(B0,L,zvec,z0).^2,mf_Jy(B0,L,zvec,z0));
  plot(hca,mf_Bx(B0,L,zvec,z0).^2,mf_Jy(B0,L,zvec,z0),'k:');
  
  hold(hca,'off');
  
  legend(hl,legs,'location','best','box','off')
end
if 1
  hca = h(isub); isub = isub + 1;  
  [N edges mid loc] = histcn([jiy(:) jey(:)],linspace(jielim(1),jielim(2),100),linspace(jielim(1),jielim(2),100));
  N(N==0) = NaN;
  pcolor(hca,mid{1:2},log10(N(:,:,ceil(end/2)))')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(counts)';
  colormap(hca,pic_colors('candy4'))
  hca.XLabel.String = 'j_{iy}';
  hca.YLabel.String = 'j_{ey}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on');
  plot(hca,jielim,jielim,'k-',-jielim,jielim,'k-')
  hold(hca,'off');
end