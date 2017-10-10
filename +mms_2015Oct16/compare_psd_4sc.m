% Compare particle distribution between the four sc
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');

ic = 1;

h = irf_plot(11);
isub = 1;

if 0 % 1 sc
  tozoom = 1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?.tlim(tint).abs},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('B_{?,DMPA}',ic),'(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z','|B|'},[1.01 0.9]);
else % sc comp
  tozoom = 1:3;
  hca = irf_panel('Bx');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dmpaB1.x.tlim(tint),dmpaB2.x.tlim(tint),dmpaB3.x.tlim(tint),dmpaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{x}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

  hca = irf_panel('By');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dmpaB1.y.tlim(tint),dmpaB2.y.tlim(tint),dmpaB3.y.tlim(tint),dmpaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{y}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

  hca = irf_panel('Bz');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dmpaB1.z.tlim(tint),dmpaB2.z.tlim(tint),dmpaB3.z.tlim(tint),dmpaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);
end

  
clim = [-1 1];
ylim = [10 1e3];
% Electron velocity in eV
for ic = 1:4
  if 1
    hca = irf_panel(irf_ssub('edistparperp ?',ic));
    c_eval('[~,hcb] = irf_spectrogram(hca,ePSDparperp?,''log'',''donotfitcolorbarlabel'');',ic)
    hold(hca,'on')  
    hold(hca,'off')
    irf_legend(hca,'(g)',[0.99 0.98],'color','k','fontsize',12)
    set(hca,'yscale','log');
    hca.CLim = clim;
    set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    hca.YLim = ylim;
    ylabel(hca,{'E_e','(eV)'});
    colormap(hca,cn.cmap('bluered3'))
    
    hhcb(isub) = hcb;
    if 1
      hcb.YTick = 0.6*[-1 1];
      hcb.YTickLabel = {'f_{perp}','f_{a/par}'};
      hcb.YLabel.String = irf_ssub('mms {?}',ic);
      hcb.YLabel.Color = mms_colors(irf_ssub('?',ic));
    end
    
    isub = isub + 1;    
  end
end
for ic = 1:4
  if 1
    hca = irf_panel(irf_ssub('edistparapar ?',ic));
    c_eval('[~,hcb] = irf_spectrogram(hca,ePSDparapar?,''log'',''donotfitcolorbarlabel'');',ic)
    irf_legend(hca,'(f)',[0.99 0.98],'color','k','fontsize',12)
    set(hca,'yscale','log');
    hca.CLim = clim;
    set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    hca.YLim = ylim;
    ylabel(hca,{'E','(eV)'});
    colormap(hca,cn.cmap('bluered3'))
    
    %hca.YLabel.String = {irf_ssub('E_{?}',ic),'(eV)'};
    
    % Colorbar labels
    hhcb(isub) = hcb;  
    if 0
      hcb.YLabel.Rotation = 0;
      hcb.YLabel.VerticalAlignment = 'middle';
      hcb.YLabel.String ={' ','f_{par}',' ','f_{apar}'};
      hcb.YLabel.Position = hcb.YLabel.Position + [1 0 0];
    end
    if 1
      hcb.YTick = 0.6*[-1 1];
      hcb.YTickLabel = {'f_{apar}','f_{par}'};
      hcb.YLabel.String = irf_ssub('mms {?}',ic);
      hcb.YLabel.Color = mms_colors(irf_ssub('?',ic));
    end
    
    isub = isub + 1;
  end
end

irf_zoom(h,'x',tint)
irf_zoom(h(tozoom),'y')
irf_plot_axis_align

%%

