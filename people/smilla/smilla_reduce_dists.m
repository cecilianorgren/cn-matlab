%% Load data
tint = irf.tint('2015-10-16T10:32:30.00Z/2015-10-16T10:34:10.00Z');
ic = 3;
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?.data./iPDistErr?.data).^2;',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('dbcsPi? = mms.get_data(''Pi_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic)
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic)

%%
% c_eval('tsXgse? = irf.ts_vec_xyz(iPDist?.time,repmat([1 0 0],iPDist?.length,1));',ic)
% c_eval('tsYgse? = irf.ts_vec_xyz(iPDist?.time,repmat([0 1 0],iPDist?.length,1));',ic)
% c_eval('tsZgse? = irf.ts_vec_xyz(iPDist?.time,repmat([0 0 1],iPDist?.length,1));',ic)
% % The skymap distributions are in their own spacecraft coordinate system,
% % so to make it properly, one needs to supply the reduce function with 
% % the wanted coordinates in the dsl coordinate system. Applied below
% c_eval('tsLdsl? = mms_dsl2gse(tsXgse?,defatt?,-1);',ic)
% c_eval('tsMdsl? = mms_dsl2gse(tsYgse?,defatt?,-1);',ic)
% c_eval('tsNdsl? = mms_dsl2gse(tsZgse?,defatt?,-1);',ic)

%% Remove some noise from ion data

c_eval('PD_counts = iPDist?_counts;',ic)
c_eval('PD_orig = iPDist?;',ic)

matEmask = PD_counts.find_low_counts('counts',6,'nMovMean',[5 5],'output','mat');
%tsElow = PD_counts.find_low_counts('counts',6,'nMovMean',[5 5],'output','energy');

PD = PD_orig.mask('energy','mat',matEmask);

%% Make movie of 2d VDFs, to get a better feeling of the turning.

%fT = 40;
%tint = time_xline_ion + 0.5*fT*[-1 1];

pdist_all = PD;
%c_eval('pdist_all_e = ePDist?;',ic)
nMean = 3;

nTimePanels = 3;
nRows = 2;
nCols = 2;
relFracTimePanels = 0.5;
[h1,h2] = initialize_combined_plot('leftright',nTimePanels,nRows,nCols,relFracTimePanels,'vertical');

%compact_panels(h,0.1,0.1)

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854

c_eval('Vi = dbcsVi?;',ic)
c_eval('scpot = scPot?;',ic)
c_eval('E = dslE?;',ic)
c_eval('B = dmpaB?;',ic)
c_eval('defatt = defatt?;',ic)

dt_scpot = scpot.time(2)-scpot.time(1);
scpot = scpot.movmean(0.03/dt_scpot); % do some moving averaging

if 1 % deflux omni i
  hca = irf_panel('deflux omni i');
  hca.ColorOrder = mms_colors('xyz');
  irf_spectrogram(hca,pdist_all.omni.deflux.specrec)
  hca.YLabel.Interpreter = 'tex';
end
if 1 % vi
  hca = irf_panel('vi');
  hca.ColorOrder = mms_colors('xyz');
  irf_plot(hca,{Vi.x,Vi.y,Vi.z},'comp')
  hca.YLabel.String = 'v_i (km/s)';
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'v_x','v_y','v_z'},[.98 0.05]);
  hca.YLabel.Interpreter = 'tex';
end

irf_zoom(h1,'x',tint)
irf_plot_axis_align(h1)

%elim = [000 Inf];
%time = time_vdf;

fontsize_leg = 9;
fontsize = 10;


vint_L = [-Inf -170];
vint_L = [-Inf Inf];
vint_M = [-Inf Inf];
vint_N = [-Inf Inf];

% Make a movie
% replace filepath with the filepath you want
directory_ = strrep(printpath,'\','');
filepath = [directory_ 'testmovie.mp4'];
vidfile = VideoWriter(filepath,'MPEG-4');
vidfile.FrameRate = 20;
open(vidfile);
     
clear F
times = pdist_all.tlim(tint).time;
for it = 10:10:times.length %1:nMovMean:times.length
  time = times(it);
  
  pd = pdist_all.tlim(time+nMean(1)*0.5*0.151*[-1 1]);
  tint_dist = [pd.time.start pd.time.stop] + 0.5*0.150*[-1 1];
  pd = pd;
  
  if exist('hmark','var'); delete(hmark); end
  c_eval('hmark = irf_pl_mark(h1,tint_dist,[0.5 0.5 0.5]);',1:numel(h1))
  
  % The coordinate system to plot in
  v1 = [1 0 0];
  v2 = [0 1 0];
  v3 = [0 0 1];
  
  % One way to get coordinates
  %c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*nMean*[-1 1]);',ic)  
  %c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*nMean*[-1 1]);',ic) 
  %c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*nMean*[-1 1]);',ic)   
  %c_eval('scaxis = mean(tsSCaxis?_lmn.tlim(pdist.time([1 end]) + 0.5*0.15*nMean*[-1 1]).data,1);',ic) 
  
  % Get those vectors in the PDist coordinate system
  t_dist_center = pd.time.start + (pd.time.stop - pd.time.start)/2; 
  tsV1 = mms_dsl2gse(irf.ts_vec_xyz(t_dist_center,v1),defatt,-1);
  tsV2 = mms_dsl2gse(irf.ts_vec_xyz(t_dist_center,v2),defatt,-1);
  tsV3 = mms_dsl2gse(irf.ts_vec_xyz(t_dist_center,v3),defatt,-1);

  V1 = tsV1.data;
  V2 = tsV2.data;
  V3 = tsV3.data;
  T = [V1; V2; V3];
 

  

  scaxis_scale = 2000;
  nSmooth = 0;
  nContours = 0;
  isub = 1;

  if 1 % f(L,M)
    hca = h2(isub); isub = isub + 1;
    vdf = pd.reduce('2D',V1,V2,'vint',vint_N,'scpot',scpot);
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    %hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
%    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    if 0 % plot B direction
      xlim = hca.XLim;
      ylim = hca.YLim;
      hold(hca,'on')
      %dt_distx = 0.030;
      %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
      B__ = B.tlim(pd.time([1 end]) + 0.5*0.03*[-1 1]);
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
    if 0 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(2)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(2)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end

    if 0 % plot eigenvectors
      hold(hca,'on')     
      eigval = tsP_eig2d_vals.resample(pd.time).data;
      EV1 = tsP_eig2d_vec1.resample(pd.time).data;
      EV2 = tsP_eig2d_vec2.resample(pd.time).data;
      quiver(hca,-EV1(1)*scaxis_scale,-EV1(2)*scaxis_scale,2*EV1(1)*scaxis_scale,2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',2)
      quiver(hca,-EV2(1)*scaxis_scale,-EV2(2)*scaxis_scale,2*EV2(1)*scaxis_scale,2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',2)
      hold(hca,'off')
      irf_legend(hca,sprintf('eig1/eig2=%.2f',eigval(1)/eigval(2)),[0.98 0.98],'color','k')
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end

  if 1 % f(L,N)
    hca = h2(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[L_vi],[N_vi]);
    vdf = pd.reduce('2D',V1,V3,'vint',vint_M,'scpot',scpot);
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 0 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  if 1 % f(M,N)
    hca = h2(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pd.reduce('2D',V2,V3,'vint',vint_L,'scpot',scpot);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 0 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(2)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(2)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
    
   

  if 1 % 1D f(vx)
    hca = h2(isub); isub = isub + 1;
    vdf = pd.reduce('1D',V1,'scpot',scpot);

    v_center = vdf.depend{1}(1,:);
    plot(hca,v_center,data1)
    hca.XLabel.String = 'v_X (km/s)';
    hca.YLabel.String = 'f_i(v_x) (s/m^4)';
    axis(hca,'square')
    
    hca.XLim = vlim*[-1 1];  
  end


  %irf_legend(h(1),{sprintf('N = %g',nMovMean)},[0.02 1.01])
  hlinks = linkprop(h2(1:3),{'CLim'});
  hlinks.Targets(1).CLim = [-10 -7.5];
  c_eval('h2(?).FontSize = 13;',1:numel(h2))
  %colormap(pic_colors('candy6'))
  colormap(irf_colormap('magma'))
  h2(1).CLim = [-10 -7.3];

  set(gcf,'color','white');

  drawnow
  pause(0.1)
  F(it) = getframe(gcf); 
  writeVideo(vidfile,F(it));  
end
close(vidfile);

