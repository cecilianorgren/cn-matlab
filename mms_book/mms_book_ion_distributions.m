mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db','/Users/cno062/Data/MMS');
%mms.db_init('local_file_db','/Volumes/mms');
db_info = datastore('mms_db');

%% Torbert event 
ic = 3;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic)
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?.data./iPDistErr?.data).^2;',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s

[enflux_new, enflux_BG, idist_new, idist_BG, Ni_new, gseVi_new, gsePi_new, ...
  Ni_bg, EnergySpectr_bg, Pres_bg, EnergySpectr_bg_self]= mms.remove_ion_penetrating_radiation_bg(iPDist3);

c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic)
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic)

% counts = (P/Perr)^2
% sum counts for all time steps,
% remove every bin which has counts<1.5

%%

L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9230];
lmn_edr = [L;M;N];

L_vi = -[-0.8906    0.4548    0.0045];
M_vi = [ 0.4539    0.8893   -0.0559];
N_vi = -[-0.0294   -0.0477   -0.9984];
lmn_vi = [L_vi; M_vi; N_vi];


L_gse = [1 0 0];
M_gse = [0 1 0];
N_gse = [0 0 1];
lmn_gse = [L_gse; M_gse; N_gse];

L_gse = [1 0 0];
M_gse = [0 1 -0.2];
%N_gse = [0 0 1];
N_gse = cross(L_gse,M_gse);
lmn_gse = [L_gse; M_gse; N_gse];

lmn = lmn_gse;
L = lmn(1,:);
M = lmn(2,:);
N = lmn(3,:);


c_eval('mvaVExB? = gseVExB?*lmn''; mvaVExB?.name = ''E LMN'';',ic)
c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn_vi''; mvaVi?.name = ''Vi LMN'';',ic)

c_eval('tsLgse? = irf.ts_vec_xyz(iPDist?.time,repmat(L,iPDist?.length,1));',ic)
c_eval('tsMgse? = irf.ts_vec_xyz(iPDist?.time,repmat(M,iPDist?.length,1));',ic)
c_eval('tsNgse? = irf.ts_vec_xyz(iPDist?.time,repmat(N,iPDist?.length,1));',ic)

c_eval('tsLdsl? = mms_dsl2gse(tsLgse?,defatt?,-1);',ic)
c_eval('tsMdsl? = mms_dsl2gse(tsMgse?,defatt?,-1);',ic)
c_eval('tsNdsl? = mms_dsl2gse(tsNgse?,defatt?,-1);',ic)

%c_eval('tsLdsl? = tsLgse?;',ic)
%c_eval('tsMdsl? = tsMgse?;',ic)
%c_eval('tsNdsl? = tsNgse?;',ic)



%% One dist at the time
%h = setup_subplots(2,2);

for dt = 15%-2:1:25
[h1,h] = initialize_combined_plot('topbottom',2,1,4,0.4,'vertical');

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854
hca = h1(isub); isub = isub + 1;
irf_plot(hca,{mvaVi3_2})
hca = h1(isub); isub = isub + 1;
irf_plot(hca,{mvaVExB3.resample(iPDist3)})
irf_zoom(h1,'x',tint_zoom)
h1(2).YLim = [-2000 2000];

isub = 1;


time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
time = time + dt;
pdist = iPDist3.tlim(time + 10*0.5*0.150*[-1 1]).elim([600 Inf]);
pdist_nobg = iPDist3_nobg.tlim(time + 7*0.5*0.150*[-1 1]).elim([300 Inf]);

c_eval('B = mvaB?.resample(pdist);',ic)  
c_eval('E = mvaE?.resample(pdist);',ic)  
c_eval('vExB = mvaVExB?.resample(pdist);',ic)  
 

c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  
nSmooth = 1;
hca = h(isub); isub = isub + 1;
vdf = pdist_nobg.reduce('2D',M_vi,L_vi);
%vdf = pdist_nobg.reduce('2D',M,L);
vdf.plot_plane(hca','smooth',nSmooth,'vectors',{dmpaB3.resample(pdist),'B'})
axis(hca,'square')
hca.XLabel.String = 'v_M (km/s)';
hca.YLabel.String = 'v_L (km/s)';
if 1 % plot B direction
  %%
  xlim = hca.XLim;
  ylim = hca.YLim;
  hold(hca,'on')
  %dt_distx = 0.030;
  %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
  B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
  B_ = mean(B__.data,1);
  B_std = std(B__.data,1);
  b = B_/norm(B_);
  B_std_inplane = std(B__.data(:,1:2),1);
  B_inplane = sqrt(sum(B_(1:2).^2));
  b = b*2000;
  quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
    B_ = B_(1:2);
end
if 0 % plot B direction
  %%
  xlim = hca.XLim;
  ylim = hca.YLim;
  hold(hca,'on')
  %dt_distx = 0.030;
  %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
  B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
  B_ = mean(B__.data,1);
  B_std = std(B__.data,1);
  b = B_/norm(B_);
  B_std_inplane = std(B__.data(:,1:2),1);
  B_inplane = sqrt(sum(B_(1:2).^2));
  if B_inplane > 2*norm(B_std_inplane)
    k = b(1)/b(2);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
    B_ = B_(1:2);
end
%    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
  
end
if 1 % plot ExB
  %%
  hold(hca,'on')
  hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.x.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ow')
  hold(hca,'off')    
end
vlim = 2500;
hca.XLim = vlim*[-1 1];
hca.YLim = vlim*[-1 1];



hca = h(isub); isub = isub + 1;
vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
%vdf = pdist_nobg.reduce('2D',[M],[N]);
vdf.plot_plane(hca','smooth',nSmooth)
axis(hca,'square')
hca.XLabel.String = 'v_M (km/s)';
hca.YLabel.String = 'v_N (km/s)';
%xlim = hca.XLim;
%ylim = hca.YLim;
if 1 % plot ExB
  %%
  hold(hca,'on')
  hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ow')
  hold(hca,'off')    
end
vlim = 2500;
hca.XLim = vlim*[-1 1];
hca.YLim = vlim*[-1 1];


hca = h(isub); isub = isub + 1;
vdf = pdist_nobg.reduce('2D',[L_vi],[N_vi]);
%vdf = pdist_nobg.reduce('2D',[L],[N]);
vdf.plot_plane(hca','smooth',nSmooth)
axis(hca,'square')
hca.XLabel.String = 'v_L (km/s)';
hca.YLabel.String = 'v_N (km/s)';
if 1 % plot ExB
  %%
  hold(hca,'on')
  hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ow')
  hold(hca,'off')    
end
vlim = 2500;
hca.XLim = vlim*[-1 1];
hca.YLim = vlim*[-1 1];



times_exact{1} = vdf.time;

if 0 % with bg
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',L,M);
  vdf.plot_plane(hca','smooth',2);
end

hca = h(isub); isub = isub + 1;
hca.ColorOrder = pic_colors('matlab');
nSmooth = 3;
iso_values = 1*10.^[-27:-15];
iso_values = [6.5e-28];
hs = pdist_nobg.plot_isosurface(hca,'vals',iso_values,'smooth',nSmooth,'fill');
%hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
%hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
axis(hca,'square')
vlim = 2000;
hca.XLim = vlim*[-1 1];
hca.YLim = vlim*[-1 1];
hca.ZLim = vlim*[-1 1];
%camlight(gca,90,-45)
%view(hca,[-1 -1 0.5])
%view(hca,[1 0.2 0.2])
view(hca,[0 1 0.2])
camlight(gca,0,0)

%hca = h(isub); isub = isub + 1;
%hs = pdist.plot_isosurface(hca,'rotate',lmn);

c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]+0.2);',1,1:numel(h1))
%cn.print(sprintf('torbert_fi_ref_dt_%g',dt))
end




%% Several at the same time.
%h = setup_subplots(2,2);
dt_all = [-15  -10 0 10 15];
dt_all = [-15:5:15]+0;
dt_all = [-15:5:15]+0;
dt_all = [-6:2:6]+0;
%dt_all = [-6:2:6]+25;
%dt_all = [-6:2:6]-00;

[h1,h] = initialize_combined_plot('topbottom',2,3,numel(dt_all),0.2,'vertical');

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854
hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
irf_plot(hca,{mvaVi3})
hca.YLabel.String = 'v_i (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
hca.YLabel.Interpreter = 'tex';

hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
irf_plot(hca,{mvaVExB3.resample(iPDist3)})
irf_zoom(h1,'x',tint_zoom)
h1(2).YLim = [-2000 2000];
hca.YLabel.String = 'v_e (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
hca.YLabel.Interpreter = 'tex';


isub = 1;
nSmooth = 1;


for dt = dt_all%(1)
%%
elim = [600 Inf];
  time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
  time = time + dt;
  tint_dist = time + 10*0.5*0.150*[-1 1];
  pdist = iPDist3.tlim(tint_dist).elim([600 Inf]);
  %pdist = iPDist3_nobg.tlim(time + 2*0.5*0.150*[-1 1]).elim([200 Inf]);
  
  counts = iPDist3_counts.tlim(tint_dist);
  counts.data(isnan(counts.data)) = 0;  
  count_sum = sum(counts.data,1);

  pdist_1crem = iPDist3.tlim(tint_dist);
  pdist_1crem.data(:,count_sum<1.5) = 0;
  pdist_1crem = pdist_1crem.elim(elim);
  pdist = pdist_1crem;

  %pdist = iPDist3_nobg.tlim(tint_dist).elim(elim);
  %pdist = iPDist3.tlim(tint_dist).elim(elim);
  if 0
    %%
    edges = -0.5:1:9;
    centers = edges(2:end)-0.5*(edges(2)-edges(1));
    N = histcounts(counts.data(:),edges);
    Nsum = histcounts(count_sum(:),edges);
    %pdist_1crem.data()
    hca = subplot(1,1,1);
    bar(centers,[N; Nsum]',2)
    legend(hca,{'Not summed','Summed'})
    hca.YScale = 'log';
  end

  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);
  
  c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
  c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)    
  
  if 0 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',[-inf inf]);
    vdf.plot_plane(hca','smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    if 1 % plot B direction
      %%
      xlim = hca.XLim;
      ylim = hca.YLim;
      hold(hca,'on')
      %dt_distx = 0.030;
      %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
      B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
      B_ = mean(B__.data,1);
      B_std = std(B__.data,1);
      b = B_/norm(B_);
      B_std_inplane = std(B__.data(:,1:2),1);
      B_inplane = sqrt(sum(B_(1:2).^2));
      b = b*2000;
      %quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
      quiver(-b(1),-b(2),2*b(1),2*b(2),0,'k','linewidth',1)
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(1:2);
    end     
    if 1 % plot ExB
      %%
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
 
  if 0 % f(L,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[L_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Ldsl],[Ndsl]);
    vdf.plot_plane(hca','smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      %%
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  if 0 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Mdsl],[Ndsl]);
    vdf.plot_plane(hca','smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      %%
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  times_exact{1} = vdf.time;
  

  nSmooth = 3;
  iso_values = 1*10.^[-27:-15];
  iso_values = [6.5e-28];
  iso_values = 7e-28;
  vlim = 3000;
  if 1 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'vals',iso_values,'smooth',nSmooth,'fill','rotate',lmn);
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')
    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[0 0 1])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.Title = [];
  end
  if 1 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'vals',iso_values,'smooth',nSmooth,'fill');
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[0 -1 0])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.Title = [];
  end
  if 1 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'vals',iso_values,'smooth',nSmooth,'fill');
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[1 0 0])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.Title = [];
  end
  
  c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]+0.2);',1,1:numel(h1))
  %cn.print(sprintf('torbert_fi_ref_dt_%g',dt))
end


c_eval('h(?).FontSize = 8;',1:numel(h))


colormap(pic_colors('candy_gray'))

hlinks_LM = linkprop(h(1:4:end),{'CLim'});
hlinks_LN = linkprop(h(2:4:end),{'CLim'});
hlinks_MN = linkprop(h(3:4:end),{'CLim'});

compact_panels(h,0.04,0.01)
drawnow
hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
delete(hb(1:end-3))

c_eval('h(?).YTickLabel = [];',4:numel(h))
c_eval('h(?).YLabel = [];',4:numel(h))
c_eval('h(?).XTick = -2000:1000:2000;',1:numel(h))
c_eval('h(?).YTick = -2000:1000:2000;',1:numel(h))
c_eval('h(?).XTickLabelRotation = 0;',1:numel(h))


c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).GridLineWidth = 1;',1:numel(h))

c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

%c_eval('h1(?).Position(3) = 0.3;',1:numel(h1))
c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

h1(1).Position = [0.1700    0.8530    0.3000    0.0848];
h1(2).Position = [0.1700    0.7672    0.3000    0.0848];


%% One time but with more overview panels
dt_all = [-30:2:30];

[h1,h] = initialize_combined_plot('leftright',3,2,2,0.4,'vertical');
%[h1,h] = initialize_combined_plot('topbottom',2,3,numel(dt_all),0.2,'vertical');

isub = 1;

tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854
hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
irf_plot(hca,{mvaVi3})
hca.YLabel.String = 'v_i (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))},[0.01 1.02]);
hca.YLabel.Interpreter = 'tex';

hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
irf_plot(hca,{mvaVExB3.resample(iPDist3)})
irf_zoom(h1,'x',tint_zoom)
h1(2).YLim = [-2000 2000];
hca.YLabel.String = 'v_e (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
hca.YLabel.Interpreter = 'tex';


 
for dt = dt_all%(1)
%%
  isub = 1;
  nSmooth = 1;
  elim = [600 Inf];
  time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
  time = time + dt;
  tint_dist = time + 10*0.5*0.150*[-1 1];
  pdist = iPDist3.tlim(tint_dist).elim([600 Inf]);
  %pdist = iPDist3_nobg.tlim(time + 2*0.5*0.150*[-1 1]).elim([200 Inf]);
  
  counts = iPDist3_counts.tlim(tint_dist);
  counts.data(isnan(counts.data)) = 0;  
  count_sum = sum(counts.data,1);

  pdist_1crem = iPDist3.tlim(tint_dist);
  pdist_1crem.data(:,count_sum<1.5) = 0;
  pdist_1crem = pdist_1crem.elim(elim);
  pdist = pdist_1crem;

  

  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',[-inf inf]);
    vdf.plot_plane(hca','smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    if 1 % plot B direction
      %%
      xlim = hca.XLim;
      ylim = hca.YLim;
      hold(hca,'on')
      %dt_distx = 0.030;
      %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
      B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
      B_ = mean(B__.data,1);
      B_std = std(B__.data,1);
      b = B_/norm(B_);
      B_std_inplane = std(B__.data(:,1:2),1);
      B_inplane = sqrt(sum(B_(1:2).^2));
      b = b*2000;
      %quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
      quiver(-b(1),-b(2),2*b(1),2*b(2),0,'k','linewidth',1)
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(1:2);
    end     
    if 1 % plot ExB
      %%
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
 
  if 1 % f(L,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[L_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Ldsl],[Ndsl]);
    vdf.plot_plane(hca','smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      %%
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Mdsl],[Ndsl]);
    vdf.plot_plane(hca','smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      %%
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  times_exact{1} = vdf.time;
  

  nSmooth = 3;
  iso_values = 1*10.^[-27:-15];
  iso_values = [6.5e-28];
  iso_values = 7e-28;
  vlim = 2500;
  if 1 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'vals',iso_values,'smooth',nSmooth,'fill','rotate',lmn);
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')
    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[0 0 1])
    view(hca,[2 -1 0.2])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.Title = [];
  end
  colormap(pic_colors('candy_gray'))
  c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]+0.2);',1,1:numel(h1))
  drawnow
  cn.print(sprintf('torbert_fi_ref_dt_%g',dt))
end