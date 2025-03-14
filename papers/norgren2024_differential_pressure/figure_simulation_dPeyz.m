% no02m = PIC('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
% ds100 = PICDist('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');
if 0 % integrate particle trajectories
  ds = ds100.twpelim(twpe).findx(103.4).findz(0.15*[-1 1]);
  for id = 1:2
    fpeaks = struct([]);
    nstep = 3;
    fxyz = ds.fxyz(1,id,[4 6],1);
    npeaks = 30;
    ff = fxyz.f;
    for ip = 1:npeaks
      [~,ind] = max(ff(:));      
      [i1,i2] = ind2sub(size(ff),ind);
      fpeaks(ip).i1 = i1;
      fpeaks(ip).i2 = i2;
      fpeaks(ip).vx = 0;
      fpeaks(ip).vy = fxyz.v(i1);
      fpeaks(ip).vz = fxyz.v(i2);
      fpeaks(ip).x = mean(fxyz.x);
      fpeaks(ip).y = 0;
      fpeaks(ip).z = mean(fxyz.z);
      ff(i1+(-nstep:1:nstep),i2+(-nstep:1:nstep)) = NaN;
    end
    fpeaks_cell{id} = fpeaks;
  end

  % integrate trajectories
  ntr_pre = 0;
  for id = 1:2
    %%
    ntr_pre = trp.ntr;
    npeaks = numel(fpeaks_cell{id});
    tspan = [16000 20000 22000]/200;
    m = 1;
    q = -1;
    for ip = 1:npeaks      
      r0 = [fpeaks(ip).x,  fpeaks(ip).y,  fpeaks(ip).z];
      v0 = [fpeaks(ip).vx, fpeaks(ip).vy, fpeaks(ip).vz];
      tr_tmp = no02m.integrate_trajectory(r0,v0,tspan,m,q);
      plot(h(1),tr_tmp.x,tr_tmp.z,'k'); pause(0.01)
      h5write_trajs(h5path_traj,tr_tmp,'id',ntr_pre+ip)
    end
    %fpeaks_cell{id} = fpeaks;
  end
end
%%
h = setup_subplots(3,3);

if 0
nrows = 2;
ncols = 4;
h(1) = subplot(nrows,ncols,[1 2]);
h(2) = subplot(nrows,ncols,ncols+[1 2]);
h(3) = subplot(nrows,ncols,3);
h(4) = subplot(nrows,ncols,3+ncols);
h(5) = subplot(nrows,ncols,4);
h(6) = subplot(nrows,ncols,4+ncols);
end

if 1 % vertical figure
  clear h
nrows = 3;
ncols = 2;
h(1) = subplot(nrows,ncols,[1 2]);
h(2) = subplot(nrows,ncols,3);
h(3) = subplot(nrows,ncols,5);
h(4) = subplot(nrows,ncols,4);
h(5) = subplot(nrows,ncols,6);
end

if 1 % vertical figure, 3 columns
  clear h
nrows = 3;
ncols = 3;
h(1) = subplot(nrows,ncols,[1 2 3]);
h(2) = subplot(nrows,ncols,4);
h(3) = subplot(nrows,ncols,7);
h(4) = subplot(nrows,ncols,5);
h(5) = subplot(nrows,ncols,8);
h(6) = subplot(nrows,ncols,6);
h(7) = subplot(nrows,ncols,9);
end

fontsize = 16;

twpe = 20000;
ds = ds100.twpelim(twpe).findx(103.4).findz(0.15*[-1 1]);


isub = 1;
if 1 % simulation frame
  hca = h(isub); isub = isub + 1;
  xlim = [101 106];
  zlim = [-1 1];
  %[ha,hb] = no02m.twpelim(twpe).xlim(xlim+[-1 1]).zlim(zlim + [-1 1]).plot_map(hca,{'peyz'},'cmap',{pic_colors('blue_red')},'clim',{0.99*[-0.005 0.005]},'A',-0.5102+[-1:0.05:1],'smooth',3);
  if 1
    ptmp = no02m.twpelim(twpe).xlim(xlim+[-1 1]).zlim(zlim + [-1 1]);
    var = ptmp.peyz;
    Ay = ptmp.A;
    pcolor(hca,ptmp.xe,ptmp.ze,smooth2(var',3))
    shading(hca,'flat')
    hold(hca,'on')
    contour(hca,ptmp.xe,ptmp.ze,Ay',-0.5102+[-1:0.05:1],'k','LineWidth',1)
    hold(hca,'off')
    hca.CLim = 0.99*[-0.005 0.005];
    colormap(hca,pic_colors('blue_red'))
    hb = colorbar(hca);    
    xlim = xlim*10;
    zlim = zlim*10;
    hca.XLabel.String = 'x/d_e';
    hca.YLabel.String = 'z/d_e';
  end
  hb.YLabel.String = 'P_{eyz}';
  hca.XLim = xlim;
  hca.YLim = zlim;
  hold(hca,'on')
  [hax,hbox] = ds.plot_boxes(hca,'de',100);
  hbox(1).LineWidth = 1; hbox(2).LineWidth = 1;
  hl = findobj(hca,'type','line');    c_eval('hl(?).LineWidth = 1;',1:numel(hl))
  hc = findobj(hca,'type','contour'); c_eval('hc(?).LineWidth = 1;',1:numel(hc))
  hca.LineWidth = 1;
  hold(hca,'off')
  hca.Title.String = '';
  hca.Layer = 'top';
end
if 0 % simulation frame
  hca = h(isub); isub = isub + 1;
  xlim = [95 115];
  zlim = [-2 2];
  no02m.twpelim(20000).xlim(xlim+[-1 1]).zlim(zlim + [-1 1]).plot_map(hca,{'-dzPezy./ne'},'cmap',{pic_colors('blue_red')},'clim',{[-0.1 0.1]},'A',0.2,'smooth',5);
  hca.XLim = xlim;
  hca.YLim = zlim;
  hold(hca,'on')
  %ds.plot_map(hca);
  hold(hca,'off')
end
if 0 % simulation frame
  hca = h(isub); isub = isub + 1;
  xlim = [95 115];
  zlim = [-2 2];
  no02m.twpelim(20000).xlim(xlim+[-1 1]).zlim(zlim + [-1 1]).plot_map(hca,{'Ey'},'cmap',{pic_colors('blue_red')},'clim',{[-0.3 0.3]},'A',0.2,'smooth',5);
  hca.XLim = xlim;
  hca.YLim = zlim;
  hold(hca,'on')
  %ds.plot_map(hca);
  hold(hca,'off')
end
%isub = isub + 1;
nsmd = 1;
if 1 % sim - vdfs
  for id = ds.nd{1}:-1:1
    hca = h(isub); isub = isub + 1;     
    %ds.update_inds({id}).plot_map(hca,[4 6],1);
    fxyz = ds.fxyz(1,id,[4 6],1);
    %pcolor(hca,fxyz.v,fxyz.v,smooth2(log10(fxyz.f),nsmd)')
    ff = smooth2(fxyz.f,nsmd);
    %ff(ff <= 0) = NaN;
    toplot = log10(ff);
    
    pcolor(hca,fxyz.v,fxyz.v,toplot')
    shading(hca,'flat')

    [VY,VZ] = ndgrid(fxyz.v,fxyz.v);
    if 1
      n = sum(sum(fxyz.f))*fxyz.dv^2;
      fvy = fxyz.f.*VY;
      fvz = fxyz.f.*VZ;
      vy = sum(fvy(:))*(fxyz.dv^2)/n;
      vz = sum(fvz(:))*(fxyz.dv^2)/n;
    end
    hold(hca,'on')
    plot(hca,vy,vz,'ko','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','auto')
    hold(hca,'off')

    colormap(hca,pic_colors('candy4'))  
    hca.XTick = [-15:5:15];
    hca.YTick = [-15:5:15]; 
    hca.XGrid = 'off';
    hca.YGrid = 'off';
    hold(hca,'on')
    plot(hca,[0 0],hca.YLim,'color',[1 1 1]*0.)
    plot(hca,hca.XLim,[0 0],'color',[1 1 1]*0.)
    hold(hca,'off')
    hca.Box = 'on';
    hca.Layer = 'top';
    hca.XLabel.String = 'v_y/v_{A}';
    hca.YLabel.String = 'v_z/v_{A}';
    %hca.Title.String = 'log_{10}f_e(v_y,v_z)';
    %hca.Title.Position = [5 10 0];
    %hca.Title.HorizontalAlignment = 'right';
    %hca.Title.VerticalAlignment = 'top';
    irf_legend(hca,'log_{10}f_e(v_y,v_z)',[0.98 0.98],'fontweight','bold','fontsize',fontsize)
    if 0 % plot fpeaks
      fpeaks = fpeaks_cell{id};
      hold(hca,'on')
      for ip = 1:numel(fpeaks)        
        plot(hca,fpeaks(ip).vy,fpeaks(ip).vz,'ok')
      end
      hold(hca,'off')
    end
    if 0 % plot trajectories within box
      if id == 1
        x1 = ds.xi1{1}(id);
        x2 = ds.xi2{1}(id);
        z1 = ds.zi1{1}(id);
        z2 = ds.zi2{1}(id);
        tr = trp.find([trp.z0]<z1);
      hold(hca,'on')
      for ip = 1:numel(fpeaks)

        
        plot(hca,fpeaks(ip).vy,fpeaks(ip).vz,'ok')
      end
      hold(hca,'off')
    end
  end  
  end
end
if 1 % sim - vdfs, only northern or southern
  for id = ds.nd{1}:-1:1
    hca = h(isub); isub = isub + 1;     
    %ds.update_inds({id}).plot_map(hca,[4 6],1);
    fxyz = ds.fxyz(1,id,[4],1);
    %pcolor(hca,fxyz.v,fxyz.v,smooth2(log10(fxyz.f),nsmd)')
    ff = smooth2(fxyz.f,nsmd);
    %ff(ff <= 0) = NaN;
    toplot = log10(ff);
    
    pcolor(hca,fxyz.v,fxyz.v,toplot')
    shading(hca,'flat')

    [VY,VZ] = ndgrid(fxyz.v,fxyz.v);
    if 1
      n = sum(sum(fxyz.f))*fxyz.dv^2;
      fvy = fxyz.f.*VY;
      fvz = fxyz.f.*VZ;
      vy = sum(fvy(:))*(fxyz.dv^2)/n;
      vz = sum(fvz(:))*(fxyz.dv^2)/n;
    end
    hold(hca,'on')
    plot(hca,vy,vz,'ko','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','auto')
    hold(hca,'off')

    colormap(hca,pic_colors('candy4'))  
    hca.XTick = [-15:5:15];
    hca.YTick = [-15:5:15]; 
    hca.XGrid = 'off';
    hca.YGrid = 'off';
    hold(hca,'on')
    plot(hca,[0 0],hca.YLim,'color',[1 1 1]*0.)
    plot(hca,hca.XLim,[0 0],'color',[1 1 1]*0.)
    hold(hca,'off')
    hca.Box = 'on';
    hca.Layer = 'top';
    hca.XLabel.String = 'v_y/v_{A}';
    hca.YLabel.String = 'v_z/v_{A}';
    %hca.Title.String = 'log_{10}f_e(v_y,v_z)';
    %hca.Title.Position = [5 10 0];
    %hca.Title.HorizontalAlignment = 'right';
    %hca.Title.VerticalAlignment = 'top';
    irf_legend(hca,'log_{10}f_e(v_y,v_z)',[0.98 0.98],'fontweight','bold','fontsize',fontsize)
    if 0 % plot fpeaks
      fpeaks = fpeaks_cell{id};
      hold(hca,'on')
      for ip = 1:numel(fpeaks)        
        plot(hca,fpeaks(ip).vy,fpeaks(ip).vz,'ok')
      end
      hold(hca,'off')
    end
    if 0 % plot trajectories within box
      if id == 1
        x1 = ds.xi1{1}(id);
        x2 = ds.xi2{1}(id);
        z1 = ds.zi1{1}(id);
        z2 = ds.zi2{1}(id);
        tr = trp.find([trp.z0]<z1);
      hold(hca,'on')
      for ip = 1:numel(fpeaks)

        
        plot(hca,fpeaks(ip).vy,fpeaks(ip).vz,'ok')
      end
      hold(hca,'off')
    end
  end  
  end
end
if 1 % sim - vdfs
  for id = ds.nd{1}:-1:1
    hca = h(isub); isub = isub + 1;     
    ds.update_inds({id}).plot_map(hca,[4 6],1,'off-diag',1);
    if 1
      fxyz = ds.fxyz(1,id,[4 6],1);
      vz = mean(mean(no02m.twpelim(twpe).xlim(fxyz.x).zlim(fxyz.z).vz([4 6])));
      vy = mean(mean(no02m.twpelim(twpe).xlim(fxyz.x).zlim(fxyz.z).vy([4 6])));
      
      [VY,VZ] = ndgrid(fxyz.v,fxyz.v);
      if 1
        n = sum(sum(fxyz.f))*fxyz.dv^2;
        fvy = fxyz.f.*VY;
        fvz = fxyz.f.*VZ;
        vy = sum(fvy(:))*(fxyz.dv^2)/n;
        vz = sum(fvz(:))*(fxyz.dv^2)/n;
      end

      dP = fxyz.f.*(VY-vy).*(VZ-vz);
      %dP(dP == 0) = NaN;
      pcolor(hca,fxyz.v,fxyz.v,smooth2(dP,nsmd)')
      shading(hca,'flat')
      hca.CLim = max(abs(hca.CLim))*[-1 1];
      hold(hca,'on')
      plot(hca,vy,vz,'ko','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','auto')
      %contour(hca,fxyz.v,fxyz.v,smooth2(log10(fxyz.f)',nsmd),'k')
      hold(hca,'off')
    end
    colormap(hca,pic_colors('blue_red'))
    hca.XTick = [-15:5:15];
    hca.YTick = [-15:5:15]; 
    hca.XGrid = 'off';
    hca.YGrid = 'off';
    hold(hca,'on')
    plot(hca,[0 0],hca.YLim,'color',[1 1 1]*0.)
    plot(hca,hca.XLim,[0 0],'color',[1 1 1]*0.)
    hca.Box = 'on';
    hca.Layer = 'top';
    hca.XLabel.String = 'v_y/v_{A}';
    hca.YLabel.String = 'v_z/v_{A}';

    %hca.Title.String = 'dP_{eyz}(v_y,v_z)';
    %hca.Title.Position = [5 10 0];
    %hca.Title.HorizontalAlignment = 'right';
    %hca.Title.VerticalAlignment = 'top';
    irf_legend(hca,'dP_{eyz}(v_y,v_z)',[0.98 0.98],'fontweight','bold','fontsize',fontsize)
  end  
end
if 0 % sim - vdfs
  for id = 1
    hca = h(2);
    ds.update_inds({id}).plot_map(hca,[4 6],1,'off-diag',1);
      axis(hca,'square')
    if 0
      fxyz = ds.fxyz(1,id,[4 6],1);
      vz = mean(mean(no02m.twpelim(twpe).xlim(fxyz.x).zlim(fxyz.z).vz([4 6])));
      vy = mean(mean(no02m.twpelim(twpe).xlim(fxyz.x).zlim(fxyz.z).vy([4 6])));
      [VY,VZ] = ndgrid(fxyz.v,fxyz.v);
      dP = fxyz.f.*(VY-vy).*(VZ-vz);
      pcolor(hca,fxyz.v,fxyz.v,smooth2(dP,1)')
      shading(hca,'flat')
      hca.CLim = max(abs(hca.CLim))*[-1 1];
      hold(hca,'on')
      plot(hca,vy,vz,'kx')
      hold(hca,'off')
    end
    colormap(hca,pic_colors('blue_red'))
    hca.XTick = [-15:5:15];
    hca.YTick = [-15:5:15]; 
    hca.XGrid = 'off';
    hca.YGrid = 'off';
    hold(hca,'on')
    plot(hca,[0 0],hca.YLim,'color',[1 1 1]*0.7)
    plot(hca,hca.XLim,[0 0],'color',[1 1 1]*0.7)
  end  
end

if ncols == 2
colormap(h(2),pic_colors('candy4'))
colormap(h(3),pic_colors('candy4'))
colormap(h(4),pic_colors('blue_red'))
colormap(h(5),pic_colors('blue_red'))
elseif ncols == 3
colormap(h(2),pic_colors('candy4'))
colormap(h(3),pic_colors('candy4'))
colormap(h(4),pic_colors('candy4'))
colormap(h(5),pic_colors('candy4'))
colormap(h(6),pic_colors('blue_red'))
colormap(h(7),pic_colors('blue_red'))

end
%linkprop(h(3:4),{'clim','xlim','ylim'});
%linkprop(h(5:6),{'clim','xlim','ylim'});

c_eval('h(?).XLim = [-14.5 5.5]; h(?).YLim = [-10 10];',2:numel(h))
c_eval('axis(h(?),''square'');',2:numel(h))
compact_panels(h(2:end),0.02,0.02)
if ncols == 2
  h(2).CLim = [-6.5 -2];
  h(3).CLim = [-6.5 -2];
  h(4).CLim = 0.008*[-1 1];
  h(5).CLim = 0.008*[-1 1];
elseif ncols == 3
  h(2).CLim = [-6.5 -2];
  h(3).CLim = [-6.5 -2];
  h(4).CLim = [-6.5 -2];
  h(5).CLim = [-6.5 -2];
  h(6).CLim = 0.008*[-1 1];
  h(7).CLim = 0.008*[-1 1];


end

c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).LineWidth = 1;',1:numel(h))

fontsize_leg = 18;
if ncols == 2
  irf_legend(h(1),{'(a)'},[-0.09 0.98],'color','k','fontsize',fontsize_leg)
  irf_legend(h(2),{'(b)'},[-0.18 0.98],'color','k','fontsize',fontsize_leg)
  irf_legend(h(3),{'(c)'},[-0.18 0.98],'color','k','fontsize',fontsize_leg)
  irf_legend(h(4),{'(d)'},[1.05 0.98],'color','k','fontsize',fontsize_leg)
  irf_legend(h(5),{'(e)'},[1.05 0.98],'color','k','fontsize',fontsize_leg)
elseif ncols == 3
  fontweight = 'bold';
  irf_legend(h(1),{'(a)'},[0.02 0.98],'color','k','fontsize',fontsize_leg,'fontweight',fontweight)
  irf_legend(h(2),{'(b)'},[0.02 0.98],'color','k','fontsize',fontsize_leg,'fontweight',fontweight)
  irf_legend(h(3),{'(c)'},[0.02 0.98],'color','k','fontsize',fontsize_leg,'fontweight',fontweight)
  irf_legend(h(4),{'(d)'},[0.02 0.98],'color','k','fontsize',fontsize_leg,'fontweight',fontweight)
  irf_legend(h(5),{'(e)'},[0.02 0.98],'color','k','fontsize',fontsize_leg,'fontweight',fontweight)
  irf_legend(h(6),{'(f)'},[0.02 0.98],'color','k','fontsize',fontsize_leg,'fontweight',fontweight)
  irf_legend(h(7),{'(g)'},[0.02 0.98],'color','k','fontsize',fontsize_leg,'fontweight',fontweight)
end  



if 0 % make lines and fonts thicker
  hb.FontWeight = 'bold';
  c_eval('h(?).FontWeight = ''bold'';',1:numel(h))
  c_eval('h(?).LineWidth = 2;',1:numel(h))
end














