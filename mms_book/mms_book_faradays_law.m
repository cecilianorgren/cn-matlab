%% Plot
B0 = 1.0;
L = 0.25;

E0 = 0.2;
Lx = 0.3;
Lz = 0.4;

%B0 = 1;
%L = 0.25;
%E0 = 0.2;
%Lx = 0.4;
%Lz = 0.4;


Bx = @(z) B0*tanh(z/L);
Bx2 = @(z) B0;
Ey = @(x,z,t) E0*exp(-x.^2/(Lx*t)^2-z.^2/(Lz*t)^2);

Ay = @(z) B0*L*log(cosh(z/L));
Ay2 = @(z) B0*z;


x = linspace(-1.5,1.5,200); dx = x(2) - x(1);
z = linspace(-1.5,1.5,100); dz = z(2) - z(1);
[X,Z] = ndgrid(x,z);

BX = Bx(Z);
EY = Ey(X,Z,t1);
EY1 = Ey(X,Z,t1);
EY2 = Ey(X,Z,t2);
BZ = Z*0;
%EY = Ey(X,Z,t2);
dA = EY;

rotEY_x = -[diff(EY,1,2)/dz, (EY(:,end)-EY(:,end-1))/dz];
rotEY_z = [diff(EY,1,1)/dz; (EY(end,:)-EY(end-1,:))/dz];
%dEYdt = (EY2-EY1)/(t2-t1);

% For quivers.
xq = linspace(-2,2,20); dxq = xq(2) - xq(1);
zq = linspace(-2,2,20); dzq = zq(2) - zq(1);
[Xq,Zq] = ndgrid(xq,zq);


BXq = Bx(Zq);
EYq = Ey(Xq,Zq,t1);
BZq = Zq*0;
%EY = Ey(X,Z,t2);

rotEYq_x = -[(EYq(:,2)-EYq(:,1))/dzq, (EYq(:,3:end)-EYq(:,1:end-2))/(2*dzq), (EYq(:,end)-EYq(:,end-1))/dzq];
rotEYq_z = [(EYq(2,:)-EYq(1,:))/dxq; (EYq(3:end,:)-EYq(1:end-2,:))/(2*dxq); (EYq(end,:)-EYq(end-1,:))/dxq];


cmap = pic_colors('blue_red'); cmap = flipdim(cmap(1:floor(size(cmap,1)/2),:),1);
cmap = flipdim(pic_colors('thermal'),1);
%cmap = flipdim(colormap('gray'),1);
colors = pic_colors('matlab');
lwidth = 1;
Bwidth = 1;
%Eqwidth = 2;
alev = -10:0.15:10;
Ecolor = 0*[.6 0.1 0.1]+0*colors(3,:).^2;


hca = subplot(1,2,1);
legs = {};

  lwidth = 3;

if 1 % Ey
  hold(hca,'on')
  hE = surf(hca,X,Z,X*0,squeeze(EY2));
  %hcont = contourf(hca,X,Z,EY2,'linewidth',1.0,'color','k');
  shading(hca,'flat')
  hca.XGrid = 'off';
  hca.YGrid = 'off';
  view([0 0 1])
  hold(hca,'off')
  legs{end+1} = 'E_y';
end

if 1 % A0  
  hold(hca,'on')
  [hl_,hl] = contour(hca,X,Z,Ay2(Z),alev,'k-','linewidth',Bwidth);
  hold(hca,'off')  
  legs{end+1} = 'B';  
end

if 0 % B0
  hold(hca,'on')
  hB2 = quiver(hca,Xq2,Zq2,BXq2,BZq2,Sq,'color',colors(1,:),'linewidth',Eqwidth);
  hold(hca,'off')
  legs{end+1} = '\Delta B';
end
if 1 % -rotX, ADD along stream linesaq 
  hold(hca,'on')
  Sq = 0.5;
  hB2 = quiver(hca,Xq,Zq,-rotEYq_x,-rotEYq_z,Sq,'color',Ecolor,'linewidth',lwidth);
  hold(hca,'off')
  legs{end+1} = '\Delta{B}=-\Delta{t}\nabla\times E';
end
if 1 % A1
  hold(hca,'on')
  hl = contour(hca,X,Z,Ay2(Z)+dA,alev,'k--','linewidth',Bwidth);  
  hold(hca,'off')
  legs{end+1} = 'B+\Delta{B}';
end


%legend(hca,legs,'fontsize',14,'location','eastoutside','box','off')
colormap(hca,cmap)
hca.CLim = [0 E0*1.3];
hca.XLim = [-1.5 1.5];
axis equal
hca.XLim = x([1 end]);
hca.YLim = z([1 end]);
hca.Box = 'on';
hca.Layer = 'top';
%axis off
hca.XTick = [];
hca.YTick = [];
hca.XLabel.String = 'L';
hca.YLabel.String = 'N';
%hca.Position = [0.1300    0.1100    0.6023    0.8150];
hca.LineWidth = 1;
hca.FontSize = 20;
%legend(hca,{'E_y','B','B-\int \nabla\times E dt','B','-\nabla\times E'},'fontsize',16,'location','eastoutside','box','off')



hca = subplot(1,2,2);
legs = {};


if 1 % Ey
  hold(hca,'on')
  %hE = surf(hca,X,Z,X*0,squeeze(EY2));
  hcont = contourf(hca,X,Z,EY2,'linewidth',1.0,'color','k');
  shading(hca,'flat')
  hca.XGrid = 'off';
  hca.YGrid = 'off';
  view([0 0 1])
  hold(hca,'off')
  legs{end+1} = 'E_y';
end

if 1 % A0  
  hold(hca,'on')
  [hl_,hl] = contour(hca,X,Z,Ay(Z),alev,'k-','linewidth',Bwidth);
  hold(hca,'off')  
  legs{end+1} = 'B';  
end

if 0 % B0
  hold(hca,'on')
  hB2 = quiver(hca,Xq2,Zq2,BXq2,BZq2,Sq,'color',colors(1,:),'linewidth',Eqwidth);
  hold(hca,'off')
  legs{end+1} = '\Delta B';
end
if 1 % -rotX, ADD along stream linesaq
  hold(hca,'on')
  Sq = 0.5;
  hB2 = quiver(hca,Xq,Zq,-rotEYq_x,-rotEYq_z,Sq,'color',Ecolor,'linewidth',lwidth);
  hold(hca,'off')
  legs{end+1} = '\Delta{B}=-\Delta{t}\nabla\times E';
end
if 1 % A1
  hold(hca,'on')
  hl = contour(hca,X,Z,Ay(Z)+dA,alev,'k--','linewidth',Bwidth);  
  hold(hca,'off')
  legs{end+1} = 'B+\Delta{B}';
end


if 0 % B-dt*rotX, ADD along stream linesaq
  hold(hca,'on')
  hB2 = quiver(hca,Xq,Zq,BXq-rotEYq_x,BZq-rotEYq_z,Sq,'color','k','linewidth',lwidth);
  hold(hca,'off')
end


%legend(hca,legs,'fontsize',14,'location','eastoutside','box','off')
colormap(hca,cmap)
hca.CLim = [0 E0*1.3];
hca.XLim = [-1.5 1.5];
axis equal
hca.XLim = x([1 end]);
hca.YLim = z([1 end]);
hca.Box = 'on';
hca.Layer = 'top';
%axis off
hca.XTick = [];
hca.YTick = [];
hca.XLabel.String = 'L';
hca.YLabel.String = 'N';
%hca.Position = [0.1300    0.1100    0.6023    0.8150];
hca.LineWidth = 1;
hca.FontSize = 20;
compact_panels(0.02,0.005)
