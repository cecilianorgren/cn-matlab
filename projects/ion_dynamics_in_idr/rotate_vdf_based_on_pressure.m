function varargout = rotate_vdf_based_on_pressure(pdist,pressure,vel)

% I made this based on Daniels rotate function
T = squeeze(pressure);
vel = squeeze(vel);


% Rotate tensor, to find the most unequal component
theta1 = 0.5*atan((2*T(2,3))./(T(3,3)-T(2,2)));
theta2 = 0.5*atan((2*T(1,2))./(T(2,2)-T(1,1)));
%theta2 = 0.5*atan((2*T(2,1))./(T(1,1)-T(2,2)));


% Rotation matrix
%R1 = [1 0 0; 0 cos(theta1) -sin(theta1); 0 sin(theta1) cos(theta1)];
%T1 = R1*(T*transpose(R1));

%R2 = [cos(theta2) -sin(theta2) 0; sin(theta2) cos(theta2) 0; 0 0 1];
R2 = [cos(theta2) -sin(theta2) 0; sin(theta2) cos(theta2) 0; 0 0 1];
T2 = R2*(T*transpose(R2));

if T2(1,1) > T2(2,2)
  theta2 = -theta2;
  R2 = [cos(theta2) -sin(theta2) 0; sin(theta2) cos(theta2) 0; 0 0 1];
  T2 = R2*(T*transpose(R2));
end

% Construct the new distributions with Ahmads shift function.
pdist_rot1 = pdist.shift(vel, 10, eye(3), 'mms');
pdist_rot2 = pdist.shift(vel, 10, R2, 'mms');

%varargout{1} = T1;
%varargout{2} = pdist_rot1;
varargout{1} = T2;
varargout{2} = pdist_rot2;

doPlot = 1;
if doPlot
  nRows = 3;
  nCols = 2;
  %h = setup_subplots(nRows,nCols,'vertical');
  [h1,h] = initialize_combined_plot('topbottom',1,3,2,0.2,'vertical');
  h1.Position(2) = h1.Position(2)-0.05;
  h1.Position(4) = h1.Position(4)*2;

  if 0
  hca = h1(1);
  colors = mms_colors('xyza');
  set(hca,'ColorOrder',colors)  
  hh = irf_plot(hca,gseVi3,'linewidth',1);
  c_eval('hh(?).Color = colors(?,:);',1:3)
  hca.YLabel.String = 'v_i (km/s)';
  hca.YLabel.Interpreter = 'tex';

  irf_zoom(h1,'x',[gseVi3.time.start gseVi3.time.stop] + 130*[1 -1])
  if exist('hmark','var'); delete(hmark); end
  hmark = irf_pl_mark(h1,pdist.time + nMovMean*0.150*0.5*[-1 1],'k');
  end
  
  isub = 1;

  vlim = 2500*[-1 1];

  vdf = pdist.reduce('2D',[1 0 0],[0 1 0]);
%  vdf_vzneg = pdist.reduce('2D',[1 0 0],[0 1 0],'vint',[-Inf 0]);
%  vdf_vzpos = pdist.reduce('2D',[1 0 0],[0 1 0],'vint',[0 Inf]);
    
  if 1 % Original distribution
    hca = h(isub); isub = isub + 1;
    vdf.plot_plane(hca);
    if 1
      hold(hca,'on')
      quiver(hca,2000*cos(-theta2),2000*sin(-theta2),'color','k','LineWidth',2)
      hold(hca,'off')
    end
    if 1
      hold(hca,'on')
      TT = T/trace(T);
      quiver(hca,0, 0, 5000*TT(1,1), 0, 'color','r','LineWidth',2)
      quiver(hca,0, 0, 0, 5000*TT(2,2), 'color','r','LineWidth',2)
      hold(hca,'off')
    end
    axis(hca,'equal')
    axis(hca,'square')
    hca.XLim = vlim;
    hca.YLim = vlim;
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
    Tstr = {sprintf('%6.3f, %6.3f, %6.3f',T(1,1),T(1,2),T(1,3));...
            sprintf('%6.3f, %6.3f, %6.3f',T(2,1),T(2,2),T(2,3));...
            sprintf('%6.3f, %6.3f, %6.3f',T(3,1),T(3,2),T(3,3))};
    irf_legend(hca,Tstr,[0.98 0.98],'color','r','fontsize',8)
  end

  if 1 % Original distribution, shifted to bulk speed frame
    hca = h(isub); isub = isub + 1;
    vdf = pdist_rot1.reduce('2D',[1 0 0],[0 1 0]);
    vdf.plot_plane(hca);
    if 1
      hold(hca,'on')
      quiver(hca,0,0,2000*cos(-theta2),2000*sin(-theta2),'color','k','LineWidth',2)
      hold(hca,'off')
    end
    if 1
      hold(hca,'on')
      TT = T/trace(T);
      quiver(hca,0, 0, 5000*TT(1,1), 0, 'color','r','LineWidth',2)
      quiver(hca,0, 0, 0, 5000*TT(2,2), 'color','r','LineWidth',2)
      hold(hca,'off')
    end
    axis(hca,'equal')
    axis(hca,'square')
    hca.XLim = vlim;
    hca.YLim = vlim;
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
    Tstr = {sprintf('%6.3f, %6.3f, %6.3f',T(1,1),T(1,2),T(1,3));...
            sprintf('%6.3f, %6.3f, %6.3f',T(2,1),T(2,2),T(2,3));...
            sprintf('%6.3f, %6.3f, %6.3f',T(3,1),T(3,2),T(3,3))};
    irf_legend(hca,Tstr,[0.98 0.98],'color','r','fontsize',8)
  end

  if 1 % Rotated distribution
    hca = h(isub); isub = isub + 1;
    vdf = pdist_rot2.reduce('2D',[1 0 0],[0 1 0]);
    vdf.plot_plane(hca);
    if 1
      hold(hca,'on')
      quiver(hca,0,0,2000*cos(-theta2+theta2),2000*sin(-theta2+theta2),'color','k','LineWidth',2)
      hold(hca,'off')
    end
    if 1
      hold(hca,'on')
      TT = T2/trace(T2);
      quiver(hca,0, 0, 5000*TT(1,1), 0, 'color','r','LineWidth',2)
      quiver(hca,0, 0, 0, 5000*TT(2,2), 'color','r','LineWidth',2)
      hold(hca,'off')
    end
    axis(hca,'equal')
    axis(hca,'square')
    hca.XLim = vlim;
    hca.YLim = vlim;
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
    Tstr = {sprintf('%6.3f, %6.3f, %6.3f',T2(1,1),T2(1,2),T2(1,3));...
            sprintf('%6.3f, %6.3f, %6.3f',T2(2,1),T2(2,2),T2(2,3));...
            sprintf('%6.3f, %6.3f, %6.3f',T2(3,1),T2(3,2),T2(3,3))};
    irf_legend(hca,Tstr,[0.98 0.98],'color','r','fontsize',8)
  end

  if 1 % Original distribution
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('1D',[0 1 0]);
    plot(hca,vdf.depend{1},vdf.data,'linewidth',1);  
    hca.XGrid = 'on'; hca.YGrid = 'on';  
    hca.XLim = vlim;
    hca.XLabel.String = 'v_y (km/s)';
    hca.YLabel.String = 'f_i (s/m^4)';
  end
  if 1 % Original distribution
    hca = h(isub); isub = isub + 1;
    vdf = pdist_rot1.reduce('1D',[0 1 0]);
    plot(hca,vdf.depend{1},vdf.data,'linewidth',1);  
    hca.XGrid = 'on'; hca.YGrid = 'on';  
    hca.XLim = vlim;
    hca.XLabel.String = 'v_y (km/s)';
    hca.YLabel.String = 'f_i (s/m^4)';
  end
  if 1 % Original distribution
    hca = h(isub); isub = isub + 1;
    vdf = pdist_rot2.reduce('1D',[0 1 0],'vint',[-Inf 0]);
    plot(hca,vdf.depend{1},vdf.data,'linewidth',1);  
    hca.XGrid = 'on'; hca.YGrid = 'on';  
    hca.XLim = vlim;
    hca.XLabel.String = 'v_y (km/s)';
    hca.YLabel.String = 'f_i (s/m^4)';
  end

  hlinks = linkprop(h,{'CLim'});
  compact_panels(h,0.04)



  tint_print = pdist.time.utc('MMHHSS_mmm');
  %cn.print(sprintf('rotating_fvxvy_by_pressure_%s',tint_print))

end