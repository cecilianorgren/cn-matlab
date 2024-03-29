% magnetic_field_lines

Bn = 0.2;
B0 = 10;
BG = 5;

Bx = @(x,y,z) x*0 + y*0 + z*B0;
By = @(x,y,z) x*0 + y*0 + z*0 + BG;
Bz = @(x,y,z) x*Bn + y*0 + z*0;

Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z) x*0 + y*0 + z*0;
Ez = @(x,y,z) x*0 + y*0 + z*0;


B = @(x,y,z) sqrt(Bx(x,y,z).^2 + By(x,y,z).^2 + Bz(x,y,z).^2);
E = @(x,y,z) sqrt(Ex(x,y,z).^2 + Ey(x,y,z).^2 + Ez(x,y,z).^2);

xmax = 1; nx = 100;
ymax = 0; ny = 1;
zmax = 0.1; nz = 20;

x = linspace(-xmax,xmax,nx);
y = linspace(-ymax,ymax,ny);
z = linspace(-zmax,zmax,nz);
[X,Y,Z] = meshgrid(x,y,z);

%% Integrate magnetic field from model to get magnetic field line
dr = 0.05;
xyz = -[0.75*xmax 0 -zmax];
zstop = zmax;
doPlot = 1;
S = 0.1;
if doPlot
  hca = subplot(1,1,1);
  quiver3(xyz(end,1),...
          xyz(end,2),...
          xyz(end,3),...
          Bx(xyz(end,1),xyz(end,2),xyz(end,3)),...
          By(xyz(end,1),xyz(end,2),xyz(end,3)),...
          Bz(xyz(end,1),xyz(end,2),xyz(end,3)),...
          S)
    hca.XLabel.String = 'X';
    hca.YLabel.String = 'Y';
    hca.ZLabel.String = 'Z';
    hold(hca,'on')
    pause(0.1)
end
while xyz(end,3)>-zstop
  dx = dr*Bx(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  dy = dr*By(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  dz = dr*Bz(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  xyz(end+1,:) = xyz(end,:) + [dx dy dz];      
  if doPlot
    quiver3(xyz(end,1),...
            xyz(end,2),...
            xyz(end,3),...
            Bx(xyz(end,1),xyz(end,2),xyz(end,3)),...
            By(xyz(end,1),xyz(end,2),xyz(end,3)),...
            Bz(xyz(end,1),xyz(end,2),xyz(end,3)),...
            S)
    pause(0.1)
  end
end
%hold(hca,'off')

xyzB = [Bx(xyz(:,1),xyz(:,2),xyz(:,3)) By(xyz(:,1),xyz(:,2),xyz(:,3)) Bz(xyz(:,1),xyz(:,2),xyz(:,3))];

%% Plot
figure(201)
nrows = 2;
ncols = 2;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

if 1 % |B|
  hca = h(isub); isub = isub + 1;
  xplot = squeeze(X);
  zplot = squeeze(Z);
  cplot = squeeze(B(X,Y,Z));
  pcolor(hca,xplot,zplot,cplot)
  shading(hca,'flat')
end
if 1 % Bx
  hca = h(isub); isub = isub + 1;
  xplot = squeeze(X);
  zplot = squeeze(Z);
  cplot = squeeze(Bx(X,Y,Z));
  pcolor(hca,xplot,zplot,cplot)
  shading(hca,'flat')
end
if 1 % Bz
  hca = h(isub); isub = isub + 1;
  xplot = squeeze(X);
  zplot = squeeze(Z);
  cplot = squeeze(Bz(X,Y,Z));
  pcolor(hca,xplot,zplot,cplot)
  shading(hca,'flat')
end
if 1 % B line
  hca = h(isub); isub = isub + 1;
  xplot = xyz(:,1);
  yplot = xyz(:,2);
  zplot = xyz(:,3);
  plot3(xplot,yplot,zplot)  
end


%%
ic = 1;
mms_2015Nov12.Bmodel

zBstart = -30*1e3;
zBstop = 30*1e3;
dr = 100;
xyz = [0 0 zBstart];
while xyz(end,3)<zBstop
  dx = dr*Bx(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  dy = dr*By(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  dz = dr*Bz(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  xyz(end+1,:) = xyz(end,:) + [dx dy dz];      
end

xyzB = [Bx(xyz(:,1),xyz(:,2),xyz(:,3)) By(xyz(:,1),xyz(:,2),xyz(:,3)) Bz(xyz(:,1),xyz(:,2),xyz(:,3))];

% Obtain electric field along the field line
xyzE = [Ex(xyz(:,1),xyz(:,2),xyz(:,3)) Ey(xyz(:,1),xyz(:,2),xyz(:,3)) Ez(xyz(:,1),xyz(:,2),xyz(:,3))];
xEint = -cumtrapz(xyz(:,1),Ex(xyz(:,1),xyz(:,2),xyz(:,3)));
yEint = -cumtrapz(xyz(:,2),Ey(xyz(:,1),xyz(:,2),xyz(:,3)));
zEint = -cumtrapz(xyz(:,3),Ez(xyz(:,1),xyz(:,2),xyz(:,3)));
xyzEint = [xEint, yEint, zEint];

% Plot 
nrows = 4;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

colors = mms_colors('matlab');
if 1 % xyz vs z
  hca = h(isub); isub = isub + 1;  
  plot(hca,xyz(:,3)*1e-3,[xyz]*1e-3) % sqrt(sum(xyz.^2,2))
  hca.YLabel.String = 'R_{@B} (km)';
  hca.XLabel.String = 'N (km)';
  legend(hca,'L','M','N')
  hca.XLim = [min(xyz(:,3)) max(xyz(:,3))]*1e-3;
end
if 1 % xyzB vs z
  hca = h(isub); isub = isub + 1;
  set(hca,'ColorOrder',colors(1:3,:))
  plot(hca,xyz(:,3)*1e-3,xyzB*1e9,'--')
  hold(hca,'on')
  set(hca,'ColorOrder',colors(1:3,:))  
  plot(hca,zObs,obsB.data)
  hold(hca,'off')
  hca.YLabel.String = 'B_{@B} (nT)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [min(xyz(:,3)) max(xyz(:,3))]*1e-3;
end
if 1 % xyzE vs z
  hca = h(isub); isub = isub + 1;
  set(hca,'ColorOrder',colors(1:3,:))
  plot(hca,xyz(:,3)*1e-3,xyzE*1e3,'--')
  hold(hca,'on')
  set(hca,'ColorOrder',colors(1:3,:))  
  plot(hca,zObs,obsE.data)
  hold(hca,'off')
  hca.YLabel.String = 'E_{@B} (mV/m)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [min(xyz(:,3)) max(xyz(:,3))]*1e-3;
end
if 1 % xyzEint vs z
  hca = h(isub); isub = isub + 1;
  plot(hca,xyz(:,3)*1e-3,xyzEint)
  hca.YLabel.String = '\phi_{@B} (V)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [min(xyz(:,3)) max(xyz(:,3))]*1e-3;
end
if 0 % xyzE vs z
  hca = h(isub); isub = isub + 1;
end

%%

if 0
  %%
  %zz = tocolumn(linspace(zBstart,zBstop,nBz)); 
%xyzB = cumsum([Bx(zz*0,zz*0,zz) By(zz*0,zz*0,zz) Bz(zz*0,zz*0,zz)],1);
hca = subplot(1,1,1);
% plot3(hca,xyz(:,1),xyz(:,2),xyz(:,3));
% hca.XLabel.String = 'X';
% hca.YLabel.String = 'Y';
% hca.ZLabel.String = 'Z';

if 0 % LMN
  plot3(hca,xyz(:,3)*1e-3,xyz(:,1)*1e-3,xyz(:,2)*1e-3);
  hca.XLabel.String = 'N'; % Z
  hca.YLabel.String = 'L'; % X
  hca.ZLabel.String = 'M'; % Y
elseif 0 % L-M-N
  plot3(hca,xyz(:,3)*1e-3,-xyz(:,2)*1e-3,xyz(:,1)*1e-3);
  hca.XLabel.String = 'N'; % Z
  hca.YLabel.String = '-M'; % -Y
  hca.ZLabel.String = 'L'; % X
elseif 1 % L-M-N
  
  qscale = 0.5;
  plot3(hca,xyz(:,3)*1e-3,-xyz(:,2)*1e-3,xyz(:,1)*1e-3);
  hold(hca,'on')
  nq = 30;
  toplot = fix(linspace(1,size(xyz,1),nq));
  quiver3(hca,xyz(toplot,3)*1e-3,-xyz(toplot,2)*1e-3,xyz(toplot,1)*1e-3,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
    -By(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale)
  
  hca.XLabel.String = 'N'; % Z
  hca.YLabel.String = '-M'; % -Y
  hca.ZLabel.String = 'L'; % X
else % L-M-N
  xyz = xyz200; 
  qscale = 0.5;
  plot3(hca,xyz(:,3)*1e-3,-xyz(:,2)*1e-3,xyz(:,1)*1e-3);
  hold(hca,'on')
  nq = 30;
  toplot = fix(linspace(1,size(xyz,1),nq));
  quiver3(hca,xyz(toplot,3)*1e-3,-xyz(toplot,2)*1e-3,xyz(toplot,1)*1e-3,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
    -By(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale)
   
  xyz = xyz1000; plot3(hca,xyz(:,3)*1e-3,-xyz(:,2)*1e-3,xyz(:,1)*1e-3);
  toplot = fix(linspace(1,size(xyz,1),nq));
  quiver3(hca,xyz(toplot,3)*1e-3,-xyz(toplot,2)*1e-3,xyz(toplot,1)*1e-3,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
    -By(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale)
  
  xyz = xyz2000; plot3(hca,xyz(:,3)*1e-3,-xyz(:,2)*1e-3,xyz(:,1)*1e-3);
  toplot = fix(linspace(1,size(xyz,1),nq));
  quiver3(hca,xyz(toplot,3)*1e-3,-xyz(toplot,2)*1e-3,xyz(toplot,1)*1e-3,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
    -By(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale)
  
  hca.XLabel.String = 'N'; % Z
  hca.YLabel.String = '-M'; % -Y
  hca.ZLabel.String = 'L'; % X
end
axis(hca,'equal')

if 1 % add quivers to see if it's ok  
  dxyz_ = 1; % 1 km = 1000 m
  xx_ = linspace(xyz(end,1),xyz(end,1),fix(abs((xyz(1,1)-xyz(end,1))/dxyz_)));
  yy_ = linspace(xyz(end,2),xyz(end,2),fix(abs((xyz(1,2)-xyz(end,2))/dxyz_)));
  zz_ = linspace(xyz(end,3),xyz(end,3),fix(abs((xyz(1,3)-xyz(end,3))/dxyz_))); 
  [XX,YY,ZZ] = meshgrid(xx_,yy_,zz_);
  quiver(hca,XX*1e-3,YY*1e-3,ZZ*1e-3,Bx(XX,YY,ZZ),By(XX,YY,ZZ),By(XX,YY,ZZ))
elseif 0
  xx_ = linspace(xyz(1,1),xyz(end,1),10);
  yy_ = linspace(xyz(1,2),xyz(end,2),10);
  zz_ = linspace(xyz(1,3),xyz(end,3),10); 
  [XX,YY,ZZ] = meshgrid(xx_,yy_,zz_);
  quiver3(hca,ZZ*1e-3,-YY*1e-3,XX*1e-3,Bz(XX,YY,ZZ)*1e9,-By(XX,YY,ZZ)*1e9,Bx(XX,YY,ZZ)*1e9)
elseif 0
  dxyz_ = 3; % 1 km = 1000 m
  xx_ = linspace(xyz(1,1),xyz(end,1),fix(abs((xyz(1,1)-xyz(end,1))/dxyz_)));
  yy_ = linspace(xyz(1,2),xyz(end,2),fix(abs((xyz(1,2)-xyz(end,2))/dxyz_)));
  zz_ = linspace(xyz(1,3),xyz(end,3),fix(abs((xyz(1,3)-xyz(end,3))/dxyz_))); 
  [XX,YY,ZZ] = meshgrid(xx_,yy_,zz_);
  quiver3(hca,ZZ*1e-3,-YY*1e-3,XX*1e-3,Bz(XX,YY,ZZ),-By(XX,YY,ZZ),Bx(XX,YY,ZZ))
end
end

