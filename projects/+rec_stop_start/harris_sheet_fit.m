%% Specify time and spacecraft
units = irf_units;
irf.log('critical')
ic = 1:4;

localuser = datastore('local','user');

%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
db_info = datastore('mms_db');   

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
if 1%all(ic == [1:4])
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
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
mf_L = @(Bx,B0,Jy) -(B0.*((Bx./B0).^2 - 1))./Jy;

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

J_edges = linspace(-2,12,70)*1e-9; % A/m2
B_edges = linspace(0,24,71)*1e-9;  % T
[BB,JJ] = ndgrid(B_edges.^2,J_edges);
B0 = 21e-9; % T
L = mf_L(BB,B0,JJ);


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
hold(hca,'on')
clim = hca.CLim;
[c,hh] = contour(hca,B_edges.^2,J_edges,L',[5000:1000:20000]);
clabel(c,hh)
hca.CLim = clim;
hold(hca,'off')

if 0
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
end

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

%% For paper, preliminary
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
mf_L = @(Bx,B0,Jy) -(B0.*((Bx./B0).^2 - 1))./Jy;

mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);
mf_JxBz = matlabFunction(JxBz);

mf_L_B = solve(Bx,L);
mf_L_J = solve(Jy,L);

%cmap_time = ;
%cmap_count = pic_colors('');
npanels = 4;
nrows = 2;
ncols = 1;

gca = gcf;
doResize = 0;
if not(isempty(gca))
  fig_position = get(gca, 'Position');
  doResize = 1;
end
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');

if doResize
  fig = h1.Parent;
  set(fig,'position',fig_position)
end
iisub = 0;
cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];

L = [1000 1500 2000 2500 4000]*1e3;
B0 = 21e-9;
J_edges = linspace(-2,12,70);
B_edges = linspace(0,22,71);
  

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{gseBav.x,gseBav.y,gseBav.z},'comp');
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
hmark = irf_pl_mark(h1(1),tint_harris);

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
  ts_yy = gseJcurl; y_str = '<J_{y}^{curl}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  %hcb.YLabel.String = 'log_{10}(counts)';
  hcb.YLabel.String = 'counts';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  
  hca_pos = hca.Position;
  hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  hlegs.Title.String = 'L = ';  
  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
end
if 1 % Bx,Jiy, time
  hca = h2(isub); isub = isub + 1;
  timeline = tint_harris(1):0.02:tint_harris(2);
  ts_yy = gseJcurl.resample(timeline); y_str = '<J_{y}^{curl}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  start_time = irf_time(ts_yy.resample(ts_yy).tlim(tint_harris).time(1),'EpochTT>utc_HH:MM:SS');
  scatter(hca,xx.^2,yy,1,tt)  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('seconds since %s',start_time);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L  
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  %hca_pos = hca.Position;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  %hlegs.Title.String = 'L = ';  
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
  
  hca.XLim = B_edges([1 end]).^2;
  hca.YLim = J_edges([1 end]);
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

h2(2).YTick = h2(1).YTick;
h2(2).XTick = h2(1).XTick;


% Clone colorbar for time and place on top of panel 1.
hca = h1(1);
h_pos = hca.Position;
hb = colorbar('peer',hca,'location','northoutside');
hb.YTick = [];


xlim = hca.XLim;
xmark = [min(hmark.XData) max(hmark.XData)];
x1_rel = xmark(1)/diff(xlim);
x2_rel = xmark(2)/diff(xlim);

hb.Position(1) = hca.Position(1) + hca.Position(3)*x1_rel;
hb.Position(3) = hca.Position(3)*(x2_rel-x1_rel);
drawnow
hb.Position(2) = hca.Position(2) + hca.Position(4);

%% For paper, preliminary 1
tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');
units = irf_units;

% Parameters for Harris fit.
syms z z0 L B0
Bx = B0*tanh((z-z0)/L);
Jy = gradient(Bx,z)/units.mu0; % -(B0*(tanh((z - z0)/L)^2 - 1))/L
JxBz = -Jy*Bx;

% Jy = -(B0*(tanh((z - z0)/L)^2 - 1))/L 
%    = -(B0*((Bx/B0)^2 - 1))/L
%     
% L  = -(B0*((Bx/B0)^2 - 1))/Jy
mf_L = @(Bx,B0,Jy) -(B0.*((Bx./B0).^2 - 1))./Jy;
mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);
mf_JxBz = matlabFunction(JxBz);
% mf_L_B = solve(Bx,L);
% mf_L_J = solve(Jy,L);


cmap_time = irf_colormap('waterfall');
%cmap_time = colormap('parula');
cmap_time = pic_colors('candy4');
cmap_count = pic_colors('candy4');
cmap_count = flipdim(colormap('bone'),1);
cmap_count = colormap('hot');
cmap_count = flipdim(pic_colors('thermal'),1);

iisub = 0;
cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];

% Parameters for current Harris current sheet
L = [1000 1500 2000 2500 4000]*1e3;
B0 = 21e-9;
J_edges = linspace(-2,12,70);
B_edges = linspace(0,22,71);
  
fontsize = 12;

% Prepare figure
npanels = 4;
nrows = 3;
ncols = 2;

h1 = irf_plot(2);
h1(1).Position = [0.15 0.75 0.7 0.2];
h1(2).Position = [0.15 0.55 0.7 0.2];
h2(1) = axes('position',[0.15 0.12 0.3 0.3]);
h2(2) = axes('position',[0.55 0.12 0.3 0.3]);

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{gseBav.x,gseBav.y,gseBav.z},'comp');
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 0 % J, Jeav, Jiav, Jav ,curl
  for comp = ['y']
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',fontsize);  
  end
end
if 1 % J, Jeav, Jiav, Jav
  for comp = ['y']
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp)},[0.02 0.99],'fontsize',fontsize);  
  end
end

irf_zoom(h1,'x',tint)
irf_zoom(h1(zoomy),'y')
irf_plot_axis_align
hmark = irf_pl_mark(h1(1),tint_harris);

% Non-TSeries panels.
isub = 1;

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
  ts_yy = gseJcurl; y_str = 'J_{y}^{curl}';
  dt = ts_yy.time(2)-ts_yy.time(1);
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},(N*dt)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  colormap(hca,cmap_count)
  %hcb.YLabel.String = 'log_{10}(counts)';
  hcb.YLabel.String = 'time (s)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  
  hca_pos = hca.Position;
%  hlegs = legend(hl,legs,'fontsize',fontsize,'location','northoutside','box','off','orientation','horizontal');
%  hlegs.Title.String = 'L = ';
  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
end
if 1 % Bx,Jiy, time
  hca = h2(isub); isub = isub + 1;
  timeline = tint_harris(1):0.05:tint_harris(2);
  ts_yy = gseJcurl.resample(timeline); 
  %ts_yy = gseJcurl;
  %y_str = '<J_{y}^{curl}>';
  y_str = 'J_{y}^{curl}';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  start_time = irf_time(ts_yy.resample(ts_yy).tlim(tint_harris).time(1),'EpochTT>utc_HH:MM:SS');
  scatter(hca,xx.^2,yy,4,tt,'o')  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('seconds since %s',start_time);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  colormap(hca,cmap_time)
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L  
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  %hca_pos = hca.Position;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  %hlegs.Title.String = 'L = ';  
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
  
  hca.XLim = B_edges([1 end]).^2;
  hca.YLim = J_edges([1 end]);
  hca.Box = 'on';
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

h2(2).YTick = h2(1).YTick;
h2(2).XTick = h2(1).XTick;


% Clone colorbar for time and place on top of panel 1.
hca = h1(1);
h_pos = hca.Position;
hb = colorbar('peer',hca,'location','northoutside');
hb.YTick = [];
colormap(hca,cmap_time)


xlim = hca.XLim;
xmark = [min(hmark.XData) max(hmark.XData)];
x1_rel = xmark(1)/diff(xlim);
x2_rel = xmark(2)/diff(xlim);

hb.Position(1) = hca.Position(1) + hca.Position(3)*x1_rel;
hb.Position(3) = hca.Position(3)*(x2_rel-x1_rel);
drawnow
hb.Position(2) = hca.Position(2) + hca.Position(4);

h1(1).FontSize = fontsize;
h1(2).FontSize = fontsize;
h2(1).FontSize = fontsize;
h2(2).FontSize = fontsize;

%hlegs.Position(1) = 0.25;

fontsize_labels = 16;

irf_legend(h1(1),'(a)',[-0.05 0.98],'fontsize',fontsize_labels);
irf_legend(h1(2),'(b)',[-0.05 0.98],'fontsize',fontsize_labels);
irf_legend(h2(1),'(c)',[-0.1 0.98],'fontsize',fontsize_labels);
irf_legend(h2(2),'(d)',[-0.1 0.98],'fontsize',fontsize_labels);

%% For paper, preliminary 1, LMN, incl viperp and VexB
tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');
units = irf_units;
ic = 1;
cs = 'lmn';

% Parameters for Harris fit.
syms z z0 L B0
Bx = B0*tanh((z-z0)/L);
Jy = gradient(Bx,z)/units.mu0; % -(B0*(tanh((z - z0)/L)^2 - 1))/L
JxBz = -Jy*Bx;

% Jy = -(B0*(tanh((z - z0)/L)^2 - 1))/L 
%    = -(B0*((Bx/B0)^2 - 1))/L
%     
% L  = -(B0*((Bx/B0)^2 - 1))/Jy
mf_L = @(Bx,B0,Jy) -(B0.*((Bx./B0).^2 - 1))./Jy;
mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);
mf_JxBz = matlabFunction(JxBz);
% mf_L_B = solve(Bx,L);
% mf_L_J = solve(Jy,L);


cmap_time = irf_colormap('waterfall');
%cmap_time = colormap('parula');
cmap_time = pic_colors('candy4');
cmap_count = pic_colors('candy4');
cmap_count = flipdim(colormap('bone'),1);
cmap_count = colormap('hot');
cmap_count = flipdim(pic_colors('thermal'),1);

iisub = 0;
cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];

% Parameters for current Harris current sheet
L = [1000 1500 2000 2500 4000]*1e3;
B0 = 21e-9;
J_edges = linspace(-2,12,70);
B_edges = linspace(0,22,71);
  
fontsize = 12;

% Prepare figure
npanels = 4;
nrows = 3;
ncols = 2;

h1 = irf_plot(5);
ymin = 0.35; ymax = 0.95; dy = (ymax-ymin)/5;

h1(1).Position = [0.15 ymin+4*dy 0.7 dy];
h1(2).Position = [0.15 ymin+3*dy 0.7 dy];
h1(3).Position = [0.15 ymin+2*dy 0.7 dy];
h1(4).Position = [0.15 ymin+1*dy 0.7 dy];
h1(5).Position = [0.15 ymin+0*dy 0.7 dy];

dy = 0.2;
h2(1) = axes('position',[0.15 0.12 0.3 dy]);
h2(2) = axes('position',[0.55 0.12 0.3 dy]);

if 1 % B lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{lmnBav.x,lmnBav.y,lmnBav.z},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end 
if 1 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end
if 1 % J, Jeav, Jiav, Jav ,curl
  comp_str = 'M';
  for comp = ['y']
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{lmnJav.(comp).resample(timeline),lmnJiav.(comp).resample(timeline),lmnJeav.(comp).resample(timeline),lmnJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp_str),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp_str),sprintf('J_{i%s}^{FPI}',comp_str),sprintf('J_{e%s}^{FPI}',comp_str),sprintf('J_{%s}^{curl}',comp_str)},[0.02 0.99],'fontsize',fontsize);  
  end
end
irf_zoom(h1,'x',tint)
irf_zoom(h1(zoomy),'y')
irf_plot_axis_align
hmark = irf_pl_mark(h1(1),tint_harris);

% Non-TSeries panels.
isub = 1;

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
if 0 % Bx,Jiy
  hca = h2(isub); isub = isub + 1;
  ts_yy = gseJcurl; y_str = 'J_{y}^{curl}';
  dt = ts_yy.time(2)-ts_yy.time(1);
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},(N*dt)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  colormap(hca,cmap_count)
  %hcb.YLabel.String = 'log_{10}(counts)';
  hcb.YLabel.String = 'time (s)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  
  hca_pos = hca.Position;
  hlegs = legend(hl,legs,'fontsize',fontsize,'location','northoutside','box','off','orientation','horizontal');
  hlegs.Title.String = 'L = ';
  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
end
if 0 % BL,JiM, time
  hca = h2(isub); isub = isub + 1;
  timeline = tint_harris(1):0.05:tint_harris(2);
  ts_yy = gseJcurl.resample(timeline); 
  %ts_yy = gseJcurl;
  %y_str = '<J_{y}^{curl}>';
  y_str = 'J_{y}^{curl}';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  start_time = irf_time(ts_yy.resample(ts_yy).tlim(tint_harris).time(1),'EpochTT>utc_HH:MM:SS');
  scatter(hca,xx.^2,yy,4,tt,'o')  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('seconds since %s',start_time);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  colormap(hca,cmap_time)
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L  
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  %hca_pos = hca.Position;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  %hlegs.Title.String = 'L = ';  
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
  
  hca.XLim = B_edges([1 end]).^2;
  hca.YLim = J_edges([1 end]);
  hca.Box = 'on';
end
if 1 % BL,JiM
  hca = h2(isub); isub = isub + 1;
  ts_yy = lmnJcurl; y_str = 'J_{M}^{curl}';
  dt = ts_yy.time(2)-ts_yy.time(1);
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = lmnBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},(N*dt)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  colormap(hca,cmap_count)
  %hcb.YLabel.String = 'log_{10}(counts)';
  hcb.YLabel.String = 'time (s)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  
  hca_pos = hca.Position;
%  hlegs = legend(hl,legs,'fontsize',fontsize,'location','northoutside','box','off','orientation','horizontal');
%  hlegs.Title.String = 'L = ';
  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
end
if 1 % BL,JiM, time
  hca = h2(isub); isub = isub + 1;
  timeline = tint_harris(1):0.05:tint_harris(2);
  ts_yy = lmnJcurl.resample(timeline); 
  %ts_yy = gseJcurl;
  %y_str = '<J_{y}^{curl}>';
  y_str = 'J_{M}^{curl}';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = lmnBav.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  start_time = irf_time(ts_yy.resample(ts_yy).tlim(tint_harris).time(1),'EpochTT>utc_HH:MM:SS');
  scatter(hca,xx.^2,yy,4,tt,'o')  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('seconds since %s',start_time);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  colormap(hca,cmap_time)
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L  
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  %hca_pos = hca.Position;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  %hlegs.Title.String = 'L = ';  
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
  
  hca.XLim = B_edges([1 end]).^2;
  hca.YLim = J_edges([1 end]);
  hca.Box = 'on';
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

h2(2).YTick = h2(1).YTick;
h2(2).XTick = h2(1).XTick;


% Clone colorbar for time and place on top of panel 1.
hca = h1(1);
h_pos = hca.Position;
hb = colorbar('peer',hca,'location','northoutside');
hb.YTick = [];
colormap(hca,cmap_time)


xlim = hca.XLim;
xmark = [min(hmark.XData) max(hmark.XData)];
x1_rel = xmark(1)/diff(xlim);
x2_rel = xmark(2)/diff(xlim);

hb.Position(1) = hca.Position(1) + hca.Position(3)*x1_rel;
hb.Position(3) = hca.Position(3)*(x2_rel-x1_rel);
drawnow
hb.Position(2) = hca.Position(2) + hca.Position(4);

h1(1).FontSize = fontsize;
h1(2).FontSize = fontsize;
h2(1).FontSize = fontsize;
h2(2).FontSize = fontsize;

%hlegs.Position(1) = 0.25;

fontsize_labels = 16;

irf_legend(h1(1),'(a)',[-0.05 0.98],'fontsize',fontsize_labels);
irf_legend(h1(2),'(b)',[-0.05 0.98],'fontsize',fontsize_labels);
irf_legend(h2(1),'(c)',[-0.1 0.98],'fontsize',fontsize_labels);
irf_legend(h2(2),'(d)',[-0.1 0.98],'fontsize',fontsize_labels);

hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;

%% Check pressure gradient length scale for Harris sheet


syms z z0 L B0 n n0
Bx = B0*tanh((z-z0)/L);
Jy = gradient(Bx,z)/units.mu0; % -(B0*(tanh((z - z0)/L)^2 - 1))/L
JxBz = -Jy*Bx;
n = n0*cosh((z-z0)/L)^(-2);
gradn = gradient(n,z);
Ln = n/gradn;

% Jy = -(B0*(tanh((z - z0)/L)^2 - 1))/L 
%    = -(B0*((Bx/B0)^2 - 1))/L
%     
% L  = -(B0*((Bx/B0)^2 - 1))/Jy
mf_L = @(Bx,B0,Jy) -(B0.*((Bx./B0).^2 - 1))./Jy;
mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);
mf_JxBz = matlabFunction(JxBz);

mf_n = matlabFunction(n);
mf_gradn = matlabFunction(gradn);
mf_Ln = matlabFunction(Ln);
% mf_L_B = solve(Bx,L);
% mf_L_J = solve(Jy,L);

L_ = 10000e3;
z_ = L_*linspace(0.5,3,100);
n0_ = 0.1*1e6;
B0_ = 21e-9;
z0_ = 0;

nrows = 4;
ncols = 2;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % B
  hca = h(isub); isub = isub + 1;
  plot(hca,z_*1e-3,mf_Bx(B0_,L_,z_,z0_)*1e-9)
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'B (nT)';
end
if 1 % n
  hca = h(isub); isub = isub + 1;
  plot(hca,z_*1e-3,mf_n(L_,n0_,z_,z0_)*1e-6)
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'n (cc)';
end
if 1 % gradn
  hca = h(isub); isub = isub + 1;
  plot(hca,z_*1e-3,mf_gradn(L_,n0_,z_,z0_)*1e-6*1e3)
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'grad n (cc/km ?)';
end
if 1 % Ln
  hca = h(isub); isub = isub + 1;
  plot(hca,z_*1e-3,mf_Ln(L_,z_,z0_)*1e-3)
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'L_n (km ?)';
end
if 1 % (Bx,1/Ln)
  hca = h(isub); isub = isub + 1;
  plot(hca,mf_Bx(B0_,L_,z_,z0_)*1e-9,1./mf_Ln(L_,z_,z0_)*1e-3)
  hca.XLabel.String = 'B_x (nT)';
  hca.YLabel.String = '1/L_n (km)';
end
if 1 % (z,Bx*Ln)
  hca = h(isub); isub = isub + 1;
  plot(hca,z_*1e-3,mf_Bx(B0_,L_,z_,z0_)*1e-9.*mf_Ln(L_,z_,z0_)*1e-3)
  hca.XLabel.String = 'z (km)';
  %hca.XLabel.String = 'B_x (nT)';
  hca.YLabel.String = 'Bx*L_n (...)';
end

%% (B,n), (B,Tiperp)
ic = 1;
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
mf_L = @(Bx,B0,Jy) -(B0.*((Bx./B0).^2 - 1))./Jy;

mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);
mf_JxBz = matlabFunction(JxBz);

mf_L_B = solve(Bx,L);
mf_L_J = solve(Jy,L);


% Make function fit of (B,n), (B,Ti)
options = optimset(); % 'TolFun',0.2
doPlot = 1;

xx = gseB1.resample(ni1).tlim(tint_harris).x.data;  

ffit = @(x,y0,a) y0*(1-a*x.^3);

% temperature
yy = Ti1perp.resample(ni1).tlim(tint_harris).data;
yy(yy<3000) = NaN;
params_T0 = [max(yy),0.0];
doPlot = 1;
cost_function = @(params) costfunction_artemyev(params,xx,yy,doPlot);
[params_T,FVAL,EXITFLAG,OUTPUT] = fminsearch(cost_function,params_T,options);

%%
% density
yy = ne1.resample(ni1).tlim(tint_harris).data;
params_n0 = [0.3,0.1];
doPlot = 1;
cost_function = @(params) costfunction_artemyev(params,xx,yy,doPlot);
[params_n,FVAL,EXITFLAG,OUTPUT] = fminsearch(cost_function,params_n0,options);

%



%%


%cmap_time = ;
%cmap_count = pic_colors('');
npanels = 3;
nrows = 2;
ncols = 1;

gca = gcf;
doResize = 0;
if not(isempty(gca))
  fig_position = get(gca, 'Position');
  doResize = 1;
end
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');

if doResize
  fig = h1.Parent;
  set(fig,'position',fig_position)
end
iisub = 0;
cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];

L = [1000 1500 2000 2500 4000]*1e3;
B0 = 21e-9;
J_edges = linspace(-2,12,70);
B_edges = linspace(0,22,71);
  
csys = 'lmn';

if 1 % B gse/lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',csys,csys,csys),ic)
  hca.YLabel.String = {sprintf('B_{%s}',csys),'(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{csys(1),csys(2),csys(3)},[0.98 0.9],'fontsize',12);
end 
if 1 % n
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % Ti
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ti');
  set(hca,'ColorOrder',mms_colors('1234'))
  
  c_eval('irf_plot(hca,{gseTi?.trace/3,Ti?perp,Ti?par},''comp'');',ic)
  hca.YLabel.String = {'T_i','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'T_{i,||}','T_{i,\perp}','T_{i,||}'}',[1.01 0.90],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
  %hca.YLim = [0 1e4];
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
hmark = irf_pl_mark(h1(1),tint_harris);

% Non-TSeries panels.
isub = 1;
if 0 % Bx,Jy
  hca = h2(isub); isub = isub + 1;
  plot(hca,abs(gseBav.tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'.')  
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
if 0 % Bx,Jiy
  hca = h2(isub); isub = isub + 1;
  ts_yy = gseJcurl; y_str = '<J_{y}^{curl}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  %hcb.YLabel.String = 'log_{10}(counts)';
  hcb.YLabel.String = 'counts';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  
  hca_pos = hca.Position;
  hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  hlegs.Title.String = 'L = ';  
  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
end
if 1 % Bx,n, time
  hca = h2(isub); isub = isub + 1;
  timeline = tint_harris(1):0.02:tint_harris(2);
  ts_yy = ne1.resample(timeline); y_str = 'n_{e}';
  yy = ts_yy.tlim(tint_harris).data;  
  xx = gseB1.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  start_time = irf_time(ts_yy.resample(ts_yy).tlim(tint_harris).time(1),'EpochTT>utc_HH:MM:SS');
  scatter(hca,xx.^1,yy,1,tt)  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('seconds since %s',start_time);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;

  xx_ = linspace(min(xx),max(xx),50);
  plot(hca,xx_,ffit(xx_,params_n(1),params_n(2)),'k-')
  
  %hca_pos = hca.Position;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  %hlegs.Title.String = 'L = ';  
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  
  hca.XLabel.String = 'B_x (nT^2)';
  hca.YLabel.String = sprintf('%s (cm^{-3})',y_str);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  irf_legend(hca,sprintf('y = %g(1 - %g B_x^2)',params_n(1),params_n(2)),[0.98 0.98]);
  
  hca.XLim(2) = 1.05*max(xx);
  %hca.YLim = J_edges([1 end]);
end
if 1 % Bx,Tiperp, time
  hca = h2(isub); isub = isub + 1;
  timeline = tint_harris(1):0.02:tint_harris(2);
  %ts_yy = gseTi1.trace.resample(timeline)/3; y_str = 'T_{i}';
  ts_yy = Ti1perp.resample(timeline); y_str = 'T_{i\perp}';
  yy = ts_yy.tlim(tint_harris).data;  
  xx = gseB1.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  start_time = irf_time(ts_yy.resample(ts_yy).tlim(tint_harris).time(1),'EpochTT>utc_HH:MM:SS');
  scatter(hca,xx.^1,yy,1,tt)  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('seconds since %s',start_time);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;

  xx_ = linspace(min(xx),max(xx),50);
  plot(hca,xx_,ffit(xx_,params_T(1),params_T(2)),'k-')
  %hca_pos = hca.Position;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  %hlegs.Title.String = 'L = ';  
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  
  hca.XLabel.String = 'B_x (nT^2)';
  hca.YLabel.String = sprintf('%s (eV)',y_str);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  irf_legend(hca,sprintf('y = %g(1 - %g B_x^2)',params_T(1),params_T(2)),[0.98 0.98]);
  
  hca.XLim(2) = 1.05*max(xx);
  %hca.XLim = B_edges([1 end]).^2;
  %hca.YLim = J_edges([1 end]);
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

%h2(2).YTick = h2(1).YTick;
h2(2).XTick = h2(1).XTick;


% Clone colorbar for time and place on top of panel 1.
hca = h1(1);
h_pos = hca.Position;
hb = colorbar('peer',hca,'location','northoutside');
hb.YTick = [];


xlim = hca.XLim;
xmark = [min(hmark.XData) max(hmark.XData)];
x1_rel = xmark(1)/diff(xlim);
x2_rel = xmark(2)/diff(xlim);

hb.Position(1) = hca.Position(1) + hca.Position(3)*x1_rel;
hb.Position(3) = hca.Position(3)*(x2_rel-x1_rel);
drawnow
hb.Position(2) = hca.Position(2) + hca.Position(4);

%% Timeline of Harris sheet thickness, taking into account potentially varying B0
tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');
units = irf_units;

% Parameters for Harris fit. Want an expression Bx(Jy)
syms z z0 L B0
Bx = B0*tanh((z-z0)/L);
Jy = gradient(Bx,z)/units.mu0; % -(B0*(tanh((z - z0)/L)^2 - 1))/L
JxBz = -Jy*Bx;

% Jy = -(B0*(tanh((z - z0)/L)^2 - 1))/L/mu0 
%    = -(B0*((Bx/B0)^2 - 1))/L/mu0
%     
% L  = -(B0*((Bx/B0)^2 - 1))/Jy/mu0
mf_L = @(Bx,B0,Jy) (B0./Jy/units.mu0).*(1-(Bx./B0).^2);

mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);
mf_JxBz = matlabFunction(JxBz);

mf_L_B = solve(Bx,L);
mf_L_J = solve(Jy,L);

% Calculate L directly from J and Bx
B0 = 22e-9;
dt_resample = 0.5;
dt_L = 1;
timeline = tint(1):dt_resample:tint(2);

dataL = mf_L(lmnBav.x.data*1e-9,B0,Jcurl.resample(lmnBav).y.data)*1e-3; % m -> km
dataL(abs(dataL)>prctile(dataL,95)) = NaN;
dataL(Jcurl.resample(lmnBav).y.data<0) = NaN;
tsL = irf.ts_scalar(lmnBav.time,dataL);

% Data is quite scattered, so make specrec of data with binning
%L_edges = linspace(0,max(dataL),20);
L_edges = linspace(0,10000*0.99,20);
t_edges = (tint(1)+-dt_L*0.5):dt_L:(tint(2)+dt_L*0.5);
[N,edges,mid,loc] = histcn([dataL,tsL.time-t_edges(1)], L_edges, t_edges-t_edges(1));
specrecL.p_label = 'counts';
specrecL.f_label = 'L (km)';
specrecL.p = N';
specrecL.t = irf_time(t_edges(1) + mid{2},'EpochTT>epoch');
specrecL.f = mid{1};
Lpeak = nan(numel(mid{2}),1);
for it = 1:numel(mid{2})
  [PKS,LOCS] = findpeaks(N(:,it),'NPeaks',1);
  if not(isempty(LOCS))
    Lpeak(it) = mid{1}(LOCS);
  end
end
tsLpeak = irf.ts_scalar(t_edges(1) + mid{2},Lpeak);

% Since the total pressure is changing, it is likely that B0 is also
% changing. So recalculate a varying B0 based on the total pressure.
% B^2/2mu0 + Pi + Pe = Ptot = B0^2/2*mu0 (RHS is asymptotical lobe field)
%  B0 = (2*mu0*Ptot)^0.5
tsPtot = irf.ts_scalar(gsePi1.time,gsePe1.resample(gsePi1).trace.data/3+PB1.resample(gsePi1).data+gsePi1.trace.data/3);
tsB0 = irf.ts_scalar(tsPtot.time,sqrt(tsPtot.data*1e-9*2*units.mu0)*1e9); tsB0.name = 'B_0'; tsB0.units = 'nT'; % nT

dataL_varB0 = mf_L(lmnBav.x.data*1e-9,tsB0.resample(lmnBav).data*1e-9,Jcurl.resample(lmnBav).y.data)*1e-3; % m -> km
dataL_varB0(abs(dataL_varB0)>prctile(dataL_varB0,95)) = NaN;
dataL_varB0(Jcurl.resample(lmnBav).y.data<0) = NaN;
tsL_varB0 = irf.ts_scalar(lmnBav.time,dataL_varB0);

% Data is quite scattered, so make specrec of data with binning
%L_edges = linspace(0,max(dataL_varB0),20);
L_edges = linspace(0,10000*0.99,20);
t_edges = (tint(1)+-dt_L*0.5):dt_L:(tint(2)+dt_L*0.5);
[N_varB0,edges,mid,loc] = histcn([dataL_varB0,tsL.time-t_edges(1)], L_edges, t_edges-t_edges(1));
specrecL_varB0.p_label = 'counts';
specrecL_varB0.f_label = 'L (km)';
specrecL_varB0.p = N_varB0';
specrecL_varB0.p(specrecL_varB0.p==0) = NaN;
specrecL_varB0.t = irf_time(t_edges(1) + mid{2},'EpochTT>epoch');
specrecL_varB0.f = mid{1};

Lpeak_varB0 = nan(numel(mid{2}),1);
for it = 1:numel(mid{2})
  [PKS,LOCS] = findpeaks(N_varB0(:,it),'NPeaks',1);
  if not(isempty(LOCS))
    Lpeak_varB0(it) = mid{1}(LOCS);
  end
end
tsLpeak_varB0 = irf.ts_scalar(t_edges(1) + mid{2},Lpeak_varB0);
%tsL = tsL.resample(timeline);

% Setup figure
fontsize = 15;
npanels = 7;
nrows = 2;
ncols = 1;

gca = gcf;
doResize = 0;
if not(isempty(gca))
  fig_position = get(gca, 'Position');
  doResize = 1;
end
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');
c_eval('h1(?).Position(1) = 0.10;',1:numel(h1))

if doResize
  fig = h1.Parent;
  set(fig,'position',fig_position)
end
iisub = 0;
cmap = colormap(pic_colors('candy4'));
doResample = 1;
isub = 0;
zoomy = [];

L = [1000 2000 3000 4000]*1e3;

J_edges = linspace(-2,12,70);
B_edges = linspace(0,22,71);
  

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{gseBav.x,gseBav.y,gseBav.z},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end 
if 1 % J, Jeav, Jiav, Jav ,curl
  for comp = ['y']      
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

if 1 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',fontsize);  
  hca(1).YLim = [0 0.499];
end
if 1 % B0
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B0');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{tsB0},'comp');
  hca.YLabel.String = {'B_0','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'const B_0'},[0.98 0.9],'fontsize',12);
end 
if 0 % L
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('L');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{tsL},'comp');
  hca.YLabel.String = {'L','(km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'const B_0'},[0.98 0.9],'fontsize',12);
end 
if 1 % L specrec
  isub = isub + 1;
  hca = irf_panel('L specrec');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,specrecL,''log'');',2)  
  set(hca,'yscale','lin');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);    
  colormap(hca,cmap)   
  hca.YLabel.Interpreter = 'tex';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  irf_legend(hca,{sprintf('B_0 = %.0f nT',B0*1e9)},[0.02 0.98],'color','k','fontsize',fontsize)
  hold(hca,'on')
  irf_plot(hca,tsLpeak,'k')
  hold(hca,'off')
  hca.YLabel.String = {'L','(km)'}; 
end
if 1 % L specrec
  isub = isub + 1;
  hca = irf_panel('L specrec var B0');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,specrecL_varB0,''log'');',2)  
  set(hca,'yscale','lin');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);    
  colormap(hca,cmap) 
  hca.YLabel.Interpreter = 'tex';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  irf_legend(hca,{sprintf('B_0 varying')},[0.02 0.98],'color','k','fontsize',fontsize)
  %colormap(hca,irf_colormap('waterfall'))
  colormap(hca,pic_colors('candy4'))
  hold(hca,'on')
  irf_plot(hca,tsLpeak_varB0,'k')
  hold(hca,'off')
  hca.YLabel.String = {'L','(km)'};
end
if 1 % compare L with const and varying B0
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('L const and varying');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{tsLpeak,tsLpeak_varB0},'comp');
  hca.YLabel.String = {'L','(km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'constant B_0','varying B_0'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('resampled to %g s (for noise reduction)',dt_L)},[0.02 0.98],'color','k')
end 

irf_zoom(h1,'x',tint)
irf_zoom(h1(zoomy),'y')
drawnow
h1(end).YLabel.Position(1) = -0.10;
drawnow
irf_plot_axis_align
hmark = irf_pl_mark(h1(1),tint_harris);
c_eval('h1(?).FontSize = fontsize;',1:numel(h1))

%% Non-TSeries panels.
isub = 1;
if 0 % Bx,Jy
  hca = h2(isub); isub = isub + 1;
  plot(hca,abs(gseBav.tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'.')  
end
clear hl
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
  ts_yy = gseJcurl; y_str = '<J_{y}^{curl}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},log10(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  %hcb.YLabel.String = 'log_{10}(counts)';
  hcb.YLabel.String = 'counts';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  
  hca_pos = hca.Position;
  hlegs = legend(hl,legs,'fontsize',fontsize,'location','northoutside','box','off','orientation','horizontal');
  hlegs.Title.String = 'L = ';  
  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
  colormap(hca,irf_colormap('waterfall'))
end
if 1 % Bx,Jiy, time
  hca = h2(isub); isub = isub + 1;
  timeline = tint_harris(1):0.02:tint_harris(2);
  ts_yy = gseJcurl.resample(timeline); y_str = '<J_{y}^{curl}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  start_time = irf_time(ts_yy.resample(ts_yy).tlim(tint_harris).time(1),'EpochTT>utc_HH:MM:SS');
  scatter(hca,xx.^2,yy,1,tt)  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('seconds since %s',start_time);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L  
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  %hca_pos = hca.Position;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  %hlegs.Title.String = 'L = ';  
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
  
  hca.XLim = B_edges([1 end]).^2;
  hca.YLim = J_edges([1 end]);
  colormap(hca,irf_colormap('waterfall'))
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

h2(2).YTick = h2(1).YTick;
h2(2).XTick = h2(1).XTick;


c_eval('h2(?).FontSize = fontsize;',1:numel(h2))

% Clone colorbar for time and place on top of panel 1.
hca = h1(1);
h_pos = hca.Position;
hb = colorbar('peer',hca,'location','northoutside');
hb.YTick = [];
colormap(hb,irf_colormap('waterfall'))


xlim = hca.XLim;
xmark = [min(hmark.XData) max(hmark.XData)];
x1_rel = xmark(1)/diff(xlim);
x2_rel = xmark(2)/diff(xlim);

hb.Position(1) = hca.Position(1) + hca.Position(3)*x1_rel;
hb.Position(3) = hca.Position(3)*(x2_rel-x1_rel);
drawnow
hb.Position(2) = hca.Position(2) + hca.Position(4);
