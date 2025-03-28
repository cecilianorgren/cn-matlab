%% Combined ion and electron distributions in xy plane
ic = 3;
time_reversal = irf_time('2017-07-11T22:34:02.641773681Z','utc>EpochTT');

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
times = times + 0.30*1;
%times = times + 0.06;
dt_dist = 2*0.061; % 0.061 is for two distributions
dt_dist = 2*0.150; % 0.061 is for two distributions

fontsize_B_amp = 13;
markersize = 5;

vint = [-Inf Inf];
%vint = [-20000 20000];
vg = -100000:2000:100000;
vg_i = -2000:50:2000;

elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;

[h1,h2] = initialize_combined_plot('topbottom',2,2,nt,0.3,'horizontal');

% Time series
if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint)*1e-3,mvaVi?.y.tlim(tint)*1e-3,mvaVi?.z.tlim(tint)*1e-3},''comp'');',ic)  
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
  c_eval('irf_plot(hca,{-1*vte?.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
  
  hca.YLabel.String = {'u_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
%  irf_legend(hca,{['u_' comps(1)],['u_' comps(2)],['u_' comps(3)],'-v_{te}'},[0.98,0.3],'fontsize',fontsize);
end
if 1 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint)*1e-3,mvaVe?.y.tlim(tint)*1e-3,mvaVe?.z.tlim(tint)*1e-3},''comp'');',ic)  
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
  c_eval('irf_plot(hca,{-1*vte?.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
  
  hca.YLabel.String = {'u_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{['u_' comps(1)],['u_' comps(2)],['u_' comps(3)],'-v_{te}'},[0.98,0.3],'fontsize',fontsize);
end

c_eval('hm(?) = irf_pl_mark(h1(?),time_reversal);',1:numel(h1))
c_eval('hm(?).Color = [0.3 0.3 0.3];',1:numel(h1))
c_eval('hm(?).LineStyle = ''--'';',1:numel(h1))
c_eval('hm(?).LineWidth = 1;',1:numel(h1))

irf_legend(h1(1),{'Tailward'},[0.1 1.02],'color','k','fontsize',15)
irf_legend(h1(1),{'Earthward'},[0.95 1.02],'color','k','fontsize',15)
irf_legend(h1(1),{'X line'},[0.56 1.02],'color','k','fontsize',15,'horizontalalignment','center')

irf_legend(h1(1),{'I'},[0.43 0.02],'color','k','fontsize',15,'horizontalalignment','center')
irf_legend(h1(1),{'II'},[0.495 0.02],'color','k','fontsize',15,'horizontalalignment','center')
irf_legend(h1(1),{'III'},[0.56 0.02],'color','k','fontsize',15,'horizontalalignment','center')
irf_legend(h1(1),{'IV'},[0.62 0.02],'color','k','fontsize',15,'horizontalalignment','center')
irf_legend(h1(1),{'V'},[0.69 0.02],'color','k','fontsize',15,'horizontalalignment','center')


% VDFs

ql = 30;
q_lw = 1.5;
dp_clim = 0.01999*[-1 1];
ds_clim = 0.01999*[-1 1];
contour_levels = 10.^([-14 -13 -12 -11 -10 -9]+0.5);

vL_Xline = -170;
leg_time = {'I','II','III','IV','V'};
times_exact = {};
isub = 1;

nMovMean = 7;
c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim);',ic)


% Electrons
for itime = 1:times.length
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('E = mvaE?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  

  c_eval('B = mvaB?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic) 
  c_eval('ve = mvaVe?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic)   


  t_dist_center = dist.time.start + (dist.time.stop-dist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data;tsMdsl?.resample(t_dist_center).data;tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);

  c_eval('B = dmpaB?.resample(dist);',ic)  
  c_eval('ve = dbcsVe?.resample(dist);',ic)    
  ve = ve*Tdsl';
  B = B*Tdsl';
  
  vdf = dist.reduce('2D',Ldsl,Mdsl,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'contour',contour_levels);
  hc_.Colorbar.YLabel.String = sprintf('log_{10} f_e(v_%s,v_%s) (s^2/m^5)',comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('log_{10} f_e(v_%s,v_%s)',comps(1),comps(2)),'(s^2/m^5)'};
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));

  
  hca.Title.String = leg_time{itime};
  
  if 0 % plot E+vxB=0
    %%
    hold(hca,'on')
    B__ = mean(B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]).data,1);
    E__ = mean(E.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]).data,1);
    vlimit_x = (E__(2)*1e-3)/(B__(3)*1e-9)*1e-6;
    vlimit_y = (E__(2)*1e-3)/(B__(1)*1e-9)*1e-6;
    plot(hca,hca.XLim,vlimit_y*[1 1],'r--')
    plot(hca,vlimit_x*[1 1],hca.YLim,'b--')    
    hold(hca,'off')
  end
  if 1 % plot B direction quiver
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_distx = 0.030;
    %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B__ = B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);
    B_ = mean(B__.data,1);
    B_std = std(B__.data,1);
    b = B_/norm(B_);
    B_std_inplane = std(B__.data(:,1:2),1);
    B_inplane = sqrt(sum(B_(1:2).^2));
    bpl = B_;

    if B_inplane > 2*norm(B_std_inplane)
      quiver(hca,-bpl(1)*ql,-bpl(2)*ql,bpl(1)*2*ql,bpl(2)*2*ql,0,'color',[0 0 0],'linewidth',q_lw)
      %plot(hca,[-bpl(1)*ql bpl(1)*ql],[-bpl(2)*ql,bpl(2)*ql],'color',[1 0 0],'linewidth',q_lw)
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(1:2);
    end
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
  end  
  if 0 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_distx = 0.030;
    %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B__ = B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);
    B_ = mean(B__.data,1);
    B_std = std(B__.data,1);
    b = B_/norm(B_);
    B_std_inplane = std(B__.data(:,1:2),1);
    B_inplane = sqrt(sum(B_(1:2).^2));
    if B_inplane > 2*norm(B_std_inplane)
      k = b(2)/b(1);
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
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
    
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hbulk = plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ok','MarkerFaceColor','w','markersize',markersize);
    hold(hca,'off')    
  end
end

% Ions
%isub = 6;
nSmooth = 0;
for itime = 1:times.length  
  c_eval('dist = iPDist?.tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
    if 1 % f(L,M)
    hca = h2(isub); isub = isub + 1;
    vdf = dist.reduce('2D',Ldsl,Mdsl);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    vdf.plot_plane(hca,'smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    %hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
    if 0 % plot B direction
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
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    end
end

  %%
isub = 6;
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);  
  % Reduce distributions
  c_eval('dist = pdist_all?.tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('E = mvaE?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  

  c_eval('B = mvaB?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic)   

  c_eval('vExB = mvaVExB?.tlim(dist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)    
  


  t_dist_center = dist.time.start + (dist.time.stop-dist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);

  c_eval('B = dmpaB?.resample(dist);',ic)  
  c_eval('ve = dbcsVe?.resample(dist);',ic)    
  ve = ve*Tdsl';
  B = B*Tdsl';
  
  vdf = dist.reduce('2D',Ldsl,Mdsl,'vint',vint,'scpot',scpot,'vg',vg_i);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'contour',contour_levels);
  hc_.Colorbar.YLabel.String = sprintf('log_{10} f_e(v_%s,v_%s) (s^2/m^5)',comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('log_{10} f_e(v_%s,v_%s)',comps(1),comps(2)),'(s^2/m^5)'};
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));

  
  hca.Title.String = leg_time{itime};
  
  if 0 % plot E+vxB=0
    %%
    hold(hca,'on')
    B__ = mean(B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]).data,1);
    E__ = mean(E.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]).data,1);
    vlimit_x = (E__(2)*1e-3)/(B__(3)*1e-9)*1e-6;
    vlimit_y = (E__(2)*1e-3)/(B__(1)*1e-9)*1e-6;
    plot(hca,hca.XLim,vlimit_y*[1 1],'r--')
    plot(hca,vlimit_x*[1 1],hca.YLim,'b--')    
    hold(hca,'off')
  end
  if 1 % plot B direction quiver
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_distx = 0.030;
    %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B__ = B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);
    B_ = mean(B__.data,1);
    B_std = std(B__.data,1);
    b = B_/norm(B_);
    B_std_inplane = std(B__.data(:,1:2),1);
    B_inplane = sqrt(sum(B_(1:2).^2));
    bpl = B_;

    if B_inplane > 2*norm(B_std_inplane)
      quiver(hca,-bpl(1)*ql,-bpl(2)*ql,bpl(1)*2*ql,bpl(2)*2*ql,0,'color',[0 0 0],'linewidth',q_lw)
      %plot(hca,[-bpl(1)*ql bpl(1)*ql],[-bpl(2)*ql,bpl(2)*ql],'color',[1 0 0],'linewidth',q_lw)
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(1:2);
    end
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
  end  
  if 0 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_distx = 0.030;
    %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B__ = B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);
    B_ = mean(B__.data,1);
    B_std = std(B__.data,1);
    b = B_/norm(B_);
    B_std_inplane = std(B__.data(:,1:2),1);
    B_inplane = sqrt(sum(B_(1:2).^2));
    if B_inplane > 2*norm(B_std_inplane)
      k = b(2)/b(1);
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
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
    
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hbulk = plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ok','MarkerFaceColor','w','markersize',markersize);
    hold(hca,'off')    
  end
end


hlinks_e = linkprop(h2(1:5),{'CLim'});
hlinks_i = linkprop(h2(6:10),{'CLim'});

colormap(pic_colors('thermal'))
%% Formatting

% for ip = 1:2
%   for it = 1:numel(times_exact)
%     axes(h1(ip))
%     tmark = [times_exact{it}(1).epochUnix times_exact{it}(end).epochUnix] + 0.5*0.03*[-1 1]; 
%     hmark = irf_pl_mark(h1(ip),tmark,[0.5 0.5 0.5]); 
%     %hmark.FaceAlpha = 0.5;
%   end
% end
c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]); hmark(?).FaceAlpha = 0.5;',1:times.length,1:2)

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

ihsub = [1:nt-1 nt+1:(2*nt-1) 2*nt+1:(3*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;

hb(5).Position(2) = h2(5).Position(2);
hb(5).Position(4) = h2(5).Position(4);

h1(1).Position(2) = 0.79;
h1(2).Position(2) = 0.79;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt) 2*nt+2:(3*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt+nt)

%h1(1).Title.String = sprintf('MMS %g',ic);
%h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
%h1.Position(3) = h1.Position(3)*0.8;

c_eval('h2(?).XLim = 0.99*[-60 60];',1:numel(h2))
c_eval('h2(?).YLim = 0.99*[-60 60];',1:numel(h2))
%c_eval('h2(?).CLim = 0.99*[-5000 5000];',6:15)
c_eval('h2(?).Layer = ''top'';',1:numel(h2))

if 0 % Change units of data from m^-3 to kg*m^-3
for ip = 6:15 
  hsurf = findobj(h2(ip),'type','surface');
  hsurf.CData = hsurf.CData*units.me;
end
end

if 1 % make TSeries of integrated dP, dS
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  isub = 2;
  tsPP = irf.ts_scalar(times,ppsum);
  hold(h1(isub),'on')
  hpp = irf_plot(h1(isub),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(isub),'off')
  
  isub = 1;
  tsSS = irf.ts_scalar(times,sssum);
  hold(h1(isub),'on')
  hpp = irf_plot(h1(isub),tsSS,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(isub),'off')
  
  h1(2).YLabel.String = sprintf('P_{e%s} (pPa)',comps);
  h1(1).YLabel.String = ['D',sprintf('_{e%s} (pPa)',comps)];
  h1(1).YLabel.Interpreter = 'tex';
  h1(2).YLabel.Interpreter = 'tex';
end

hbb = findobj(gcf,'type','colorbar'); hbb = hbb(end:-1:1);

%hbb(2).YLabel.String = 'dP^e_{LM}/dv_Ldv_M (kg/m^3)';
%hbb(3).YLabel.String = 'dS^e_{LM}/dv_Ldv_M (kg/m^3)';

c_eval('h1(?).XGrid = ''off''; h1(?).YGrid = ''off'';',1:2)

c_eval('h2(?).Color = 1*[1 1 1];',1:15)
c_eval('h2(?).XGrid = ''off''; h2(?).YGrid = ''off'';',1:15)
color_grid = 0.5*[1 1 1];
c_eval('hold(h2(?),"on"); plot(h2(?),vg*0,vg,''color'',color_grid); plot(h2(?),vg,vg*0,''color'',color_grid);',1:15)


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};
nInd = 1;
for ii = 1:2
  irf_legend(h1(ii),legends{nInd},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  nInd = nInd + 1;
  h1(ii).FontSize = 14;
end

for ii = 1:numel(h2)
  irf_legend(h2(ii),legends{nInd},[0.04 0.97],'color',[0 0 0],'fontsize',14)
  nInd = nInd + 1;
  h2(ii).FontSize = 14;
end

