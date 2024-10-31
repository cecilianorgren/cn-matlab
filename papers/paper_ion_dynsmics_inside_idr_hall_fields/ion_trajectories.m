%no02m = PIC('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

%% Integrate trajectories
tstart = 120;
tstop = 140;
m = 1;
q = 1;
pic = no02m.twcilim(tstart);
r0 = [110 0 5];
v0 = [0 0 -0.2];
vx0 = mean(mean(pic.xlim(r0(1)+0.02*[-1 1]).zlim(r0(3)+0.02*[-1 1]).vix));
vy0 = mean(mean(pic.xlim(r0(1)+0.02*[-1 1]).zlim(r0(3)+0.02*[-1 1]).viy));
vz0 = mean(mean(pic.xlim(r0(1)+0.02*[-1 1]).zlim(r0(3)+0.02*[-1 1]).viz));
v0 = [vx0, vy0, vz0];
ptmp = pic.integrate_trajectory_constant_EB(r0,v0,tstart,tstop,m,q);


%% Figure
 
nRows = 1;
nCols = 2;
ip = 0;
for iRow = 1:nRows
  for iCol = 1:nCols
    ip = ip + 1;
    h(ip) = subplot(nRows,nCols,ip);
  end
end
isub = 1;

colors = pic_colors('matlab');
linewidth = 1;
fontsize = 13;

if 1
  hca = h(isub); isub = isub + 1;
  %x = 1:10;
  %y = 1:10;
  %z = x.^2;
  x = ptmp.x;
  y = ptmp.y;
  z = ptmp.z;
  
  plot3(hca,x,y,z,'linewidth',linewidth,'color',colors(1,:))


  hca.XLabel.String = 'L';
  hca.YLabel.String = 'M';
  hca.ZLabel.String = 'N';

  if 1 % Plot projections on the 3 planes
    hold(hca,'on')
    xlim = hca.XLim;
    ylim = hca.YLim;
    zlim = hca.ZLim;
    
    

    % x, y, z plane
    plot3(hca,x*0+xlim(2),y,z,'Color',colors(1,:).^0.5)
    plot3(hca,x,y*0+ylim(2),z,'Color',colors(1,:).^0.5)
    plot3(hca,x,y,z*0+zlim(1),'Color',colors(1,:).^0.5)
      


    hold(hca,'off')
  end
  hca.Box = 'on';
end









