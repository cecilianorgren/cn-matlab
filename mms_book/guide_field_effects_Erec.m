%% Define magnetic field model
B0 = 20e-9;
Bg = 0.5*B0;
Bn = 2e-9;
ER = 1e-3;
L = 3e3;

Bx = @(x,y,z) B0*tanh(z/L); % L
By = @(x,y,z) Bg + 0.7*B0*2*(z/L).*exp(-z.^2/L^2); % M
Bz = @(x,y,z) z*0 + Bn*x/L; % N

Ex = @(x,y,z) 0;
Ey = @(x,y,z) Er;
Ez = @(x,y,z) 0;

Babs = @(x,y,z) sqrt(Bx(x,y,z).^2 + By(x,y,z).^2 + Bz(x,y,z).^2);


zv = 5*L*linspace(-1,1,100);
xv = 5*L*linspace(-1,1,101);
xp = L;
[XV,ZV] = ndgrid(xv,zv);

h = setup_subplots(4,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,zv/L,Bx(L,0,zv)/B0,zv/L,By(L,0,zv)/B0,zv/L,Bz(L,0,zv)/B0,zv/L,Babs(L,0,zv)/B0)
irf_legend(hca,{'B_x','B_y','B_z'}',[1.01,0.98])

if 0 % E
  hca = h(isub); isub = isub + 1;
  plot(hca,zv,Ex(L,0,zv),zv,Ey(L,0,zv),zv,Ez(L,0,zv))
end

hca = h(isub); isub = isub + 1;
pcolor(hca,XV,ZV,Bx(XV,0,ZV))
shading(hca,'flat')
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'B_x';

hca = h(isub); isub = isub + 1;
pcolor(hca,XV,ZV,Bz(XV,0,ZV))
shading(hca,'flat')
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'B_z';

%% Make streamlines

zv = 5*L*linspace(-1,1,100);
xv = 5*L*linspace(-1,1,101);
yv = 0;
xp = L;
[XV,ZV] = meshgrid(xv,zv);
U = Bx(XV,0,ZV); 
V = By(XV,0,ZV); 
W = Bz(XV,0,ZV);
STARTX = -5;
STARTY = 0;
STARTZ = -2;
hs = stream2(XV,ZV,U,W,STARTX,STARTZ);

%% Integrate magnetic field from model to get magnetic field line

zBstart = -L;
zBstop = L;
dr = 100;
xyz = [-5*L 0 zBstart];
while xyz(end,3)<zBstop
  dx = dr*Bx(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  dy = dr*By(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  dz = dr*Bz(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  xyz(end+1,:) = xyz(end,:) + [dx dy dz];      
end

xyzB = [Bx(xyz(:,1),xyz(:,2),xyz(:,3)) By(xyz(:,1),xyz(:,2),xyz(:,3)) Bz(xyz(:,1),xyz(:,2),xyz(:,3))];

%% Plot 3D/2D figure

h = setup_subplots(1,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot3(hca,xyz(:,1),xyz(:,2),xyz(:,3))

