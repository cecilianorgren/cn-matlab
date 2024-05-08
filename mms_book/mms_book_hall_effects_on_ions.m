%% Define magnetic and electric field
%ic = 1; mms_2015Nov12.Bmodel;

lscale = 1; % make everything wider

a = 11e3*lscale; % thickness of current sheet, m
b = 8e3*lscale; % bifurcation length scale, m
c = 8e3*lscale;
d = 13e3*lscale;
g = 12e3*lscale;

zoffsy = 1*1e3*lscale;


limN = 30e3*1e-4*lscale;
Er = -1.0*1e-3; % reconnection electric field, V/m
Ey_inflow = 0*-1e-3;
E0 = 1*2*-1e-3;
Enas = -2.5e-3;
offsEnas = 15e3*lscale;

dzEr = 5e3;
dE = 12e3*lscale; % thickness of current sheet, m

dzE = 6e3*lscale;
BH = 4e-9;
BLoffset = 0e-9;
EL = 0*-0.5*1e-3;        


Bg = 10*1e-9; % guide field, T
Bn = 1*1e-9; % normal field, T
B0 = 5*1e-9;
Lz = 5*1e3; % m


Bx = @(x,y,z) x*0 + y*0 + z*0 - B0*tanh(z/Lz); % L
By = @(x,y,z) x*0 + y*0 + z*0 + Bg ; % M
Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn; % N
B = @(x,y,z) sqrt(Bx(x,y,z).^2 + By(x,y,z).^2 + Bz(x,y,z).^2);


Ex = @(x,y,z) x*0 + y*0 + z*0 + 0*EL*(-exp(-(z.^2)/d^2));
Ey = @(x,y,z) x*0 + y*0 + z*0 + 1*(Er*(exp(-((z-dzEr).^2)/g^2))) + Ey_inflow;
Ez = @(x,y,z) x*0 + y*0 + z*0 + 0*(E0*sin((z-dzE)/dE).*exp(-(z-dzE).^2/(2*b^2)) + 1*Enas*exp(-(z-offsEnas).^2./d.^2) + -0.8*Enas*exp(-(z+offsEnas).^2./2/d.^2));
Ez = @(x,y,z) x*0 + y*0 + z*0 + 0*(2e-3*exp(-(z+0*1e3).^2/((7*1e3)^2)) - 2*1e-3*exp(-(z-10*1e3).^2/((7*1e3)^2))); % N




h = setup_subplots(1,1,1);
isub = 1;
hca = h(isub); isub = isub + 1;

Bcolor = [0 0 0];

xyz0_all = {[0 0 zBstart]};

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
  
  nq = 20;
  qscale  = 0.2;

  toplot = fix(linspace(1,size(xyz,1),nq));
  toplot = [1 fix(size(xyz,1)*0.85)];

  hold(hca,'on')
  dL = 0;
  dM = 75; % km
  dN = 0;
  hB = plot3(hca,xyz(:,1)*1e-3 - dL,xyz(:,2)*1e-3 - dM,xyz(:,3)*1e-3 - dN,'color',Bcolor,'LineWidth',2);      
  quiver3(hca,...
     xyz(toplot,1)*1e-3 - dN,...
     xyz(toplot,2)*1e-3 - dM,...
     xyz(toplot,3)*1e-3 - dL,...
     Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     qscale,'color',Bcolor)

  dL = 0;
  dM = 10; % km
  dN = 0;
  plot3(hca,xyz(:,1)*1e-3 - dL,xyz(:,2)*1e-3 - dM,xyz(:,3)*1e-3 - dN,'color',Bcolor,'LineWidth',2);      
  quiver3(hca,...
     xyz(toplot,1)*1e-3 - dN,...
     xyz(toplot,2)*1e-3 - dM,...
     xyz(toplot,3)*1e-3 - dL,...
     Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     qscale,'color',Bcolor)
   hold(hca,'off')

%% With Ay
a = 5;
b = 2;
x = a*linspace(-10,0,1000);
y = b*linspace(-10,10,1001);
z = linspace(-10,10,5);
[X,Y] = meshgrid(x,y);
dx = x(2) - x(1);
dy = y(2) - y(1);
%dz = z(2) - z(1);
x_xline = x;
y_xline = x*b/a;

colors = pic_colors('matlab');
colors = [colors; colors(end:-1:1,:)];
linecolor_top = colors(1,:);
linecolor_bot = colors(2,:);
color_background = [0 0 0];
color_separatrix = [0 0 0];
linewidth = 3;
fontsize = 16;


Ay = @(x,y) (x/a).^2 - (y/b).^2;
AY0 = Ay(X,Y);

Astep = 20;
AY = AY0;
AYlev0 = -100:Astep:(100 + Astep); 
AYlev0 = -100:Astep:-0.5;
AYlev0 = 2:Astep:100;

hca = subplot(1,1,1);
hca.NextPlot = 'add';

S = contourcs(x,y,AY,AYlev0);
for is = 1:numel(S)
  sx = interp1(1:numel(S(is).X),S(is).X,1:0.01:numel(S(is).X));
  sy = interp1(1:numel(S(is).Y),S(is).Y,1:0.01:numel(S(is).Y));
  zmin = 0; zmax = 50;
  dsz = 0.0002;
  sz = cumsum(dsz*sign(sy).*sx);
  sz(sx<zmin) = sz(sx<zmin) + zmax;
  
  plx = sx(sy>=0); 
  ply = sy(sy>=0);
  plz = sz(sy>=0);
  plot3(hca,plx,ply,plz,'color',linecolor_top,'linewidth',linewidth)
  %plot3(hca,plx,ply,ply*0+zmax,'color',linecolor_top,'linewidth',linewidth)
  %plot3(hca,plx,ply,ply*0+zmin,'color',linecolor_top,'linewidth',linewidth)

  plx = sx(sy<=0); 
  ply = sy(sy<=0);
  plz = sz(sy<=0);
  plot3(hca,plx,ply,plz,'color',linecolor_bot,'linewidth',linewidth)
  %plot3(hca,plx,ply,plx*0+zmax,'color',linecolor_bot,'linewidth',linewidth)
  %plot3(hca,plx,ply,plx*0+zmin,'color',linecolor_bot,'linewidth',linewidth)


end

%arrow3([0 0 zmin],[0 0 zmax],'-k2')
hca.Box = 'on';
hca.BoxStyle = 'full';
hca.LineWidth = 1;
hca.FontSize = fontsize;
hca.NextPlot = 'replacechildren';
hca.View = [-13.8000   26.3072];
%axis equal

hca.XLabel.String = 'L';
hca.YLabel.String = 'N';
hca.ZLabel.String = 'M';

hca.XTick = [];
hca.YTick = [];
hca.ZTick = [];

%% With PIC
%pic = PIC('/Users/cecilia/Data/PIC/rec_onset_4/data_h5/fields_F10_E01.h5');

it = 44;%pic.nt;

%% Integrate trajectories
vy = pictmp.viy;
vx = pictmp.vix;
vz = pictmp.viz;
F_vx = griddedInterpolant(X,Z,vx);
F_vy = griddedInterpolant(X,Z,vy);
F_vz = griddedInterpolant(X,Z,vz);


x0 = [20 21 22 23 24];
z0 = [3 3 3 3 3];

ip = 1;

v0 = [0 0.00 -0.05];
vx0 = F_vx(x0',z0');
vy0 = F_vy(x0',z0');
vz0 = F_vx(x0',z0');

tstart = 0;
tstop = 150;
m = 100;
q = 1;
%clear particles
for ipart = 4:numel(x0)
  r0 = [x0(ipart) 0 z0(ipart)];
  v0 = [vx0(ipart) vy0(ipart) vz0(ipart)];
  out = pic(it).integrate_trajectory_constant_EB(r0,v0,tstart,tstop,m,q);
  particles(ipart) = out;
end
%plot3(gca,out.x,out.y,out.z)
%%
pictmp = pic(it).xlim(mean(pic.xi)+[-10 10]);
pictmp = pic(it);
Ay = pictmp.A;
By = pictmp.By;
Bx = pictmp.Bx;
Bz = pictmp.Bz;
x = pictmp.xi;
z = pictmp.zi;
[X,Z] = ndgrid(x,z);

AYlev = 10.5:.05:11.34;
AYlev = 10.5:.05:11.09;
S = contourcs(x,z,Ay',AYlev);

Fx = griddedInterpolant(X,Z,Bx);
Fy = griddedInterpolant(X,Z,By);
Fz = griddedInterpolant(X,Z,Bz);

clear B
for iLine = 1:numel(S)
  dx = diff(S(iLine).X);
  dz = diff(S(iLine).Y);

  bx = Fx(S(iLine).X,S(iLine).Y);
  by = Fy(S(iLine).X,S(iLine).Y);
  bz = Fz(S(iLine).X,S(iLine).Y);
  
  dy = dx.*by(2:end)./bx(2:end);
  y = cumsum(dy);

  %dbx = diff(bx);
  %dby = diff(by);
  %dbz = diff(bz);
  
  B(iLine).bx = bx;
  B(iLine).by = by;
  B(iLine).bz = bz;
  

  B(iLine).x = S(iLine).X(2:end);
  B(iLine).y = y;
  B(iLine).z = S(iLine).Y(2:end);
  
end




h = setup_subplots(1,1);
isub = 1;

if 0
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x,z,Fy(X,Z)')
  shading(hca,'flat')
  hcb = colorbar(hca);
  hold(hca,'on')
  for iLine = 1:numel(S)
    plot(hca,S(iLine).X,S(iLine).Y,'k')  
  end
  hold(hca,'off')
end

if 0
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x,z,By')
  shading(hca,'flat')
  hold(hca,'on')
end

hca = h(isub); isub = isub + 1;
for iLine = 1:numel(S)
  plx = B(iLine).x;
  ply = B(iLine).y;
  plz = B(iLine).z;
  irem = find(plz<-10);
  plx(irem) = [];
  ply(irem) = [];
  plz(irem) = [];
  plot3(hca,plx,ply,plz,'linewidth',1,'color',[0 0 0])
  hold(hca,'on')
end
hold(hca,'off')
if 0 % add particles
  hold(hca,'on')
  for ip = 1:numel(particles)
    part = particles(ip);
    %plot3(hca,out.x,out.y,out.z,'k')
    plot3(hca,part.x,part.y,part.z,'k')
    hold(hca,'on')
  end
  hold(hca,'off')
end

hca.XLabel.String = 'L';
hca.YLabel.String = 'M';
hca.ZLabel.String = 'N';
