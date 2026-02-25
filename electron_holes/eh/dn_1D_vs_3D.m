units = irf_units;
dn = @(z,lz,r,lr) units.eps0/units.e*(r.^2./lr.^4-2./lr.^2-1./lz.^2+z.^2./lz.^4).*exp(-z.^2/2./lz - r.^2/2./lr);

dn_xyz = @(x,lx,y,ly,z,lz) units.eps0/units.e*(3-x.^2./lx.^2-y.^2./ly.^2-z.^2./lz.^2).*exp(-x.^2/2./lx-y.^2/2./ly-z.^2/2./lz);

lz = 1;
ly = 1;
lx = 1;
lr = 1*lz;
z = linspace(-10*lz,10*lz,101);
x = z;
y = z;
r = [0 0.5 1 1.5 2]*lr;
%r = z(z>=0);
%r = r(1:3:end);
[Z,R] = meshgrid(z,r);
DN = dn(Z,lz,R,lr);
[X_xyz,Y_xyz,Z_xyz] = meshgrid(x,y,z);
DN_xyz = dn_xyz(X_xyz,lx,Y_xyz,ly,Z_xyz,lz);

DN_xz = sum(DN_xyz,2);
DN_z = sum(DN_xz,1);

nrows = 2;
ncols = 3;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;


hca = h(isub); isub = isub + 1;
plotyy(hca,z,dn(z,lz,0,100*lr),z,cumsum(dn(z,lz,0,100*lr)))

hca = h(isub); isub = isub + 1;
plotyy(hca,z,squeeze(DN_z),z,cumsum(squeeze(DN_z)))
hca.Title.String = sprintf('sum = %g',sum(squeeze(DN_z)));

hca = h(isub); isub = isub + 1;
pcolor(hca,R,Z,DN)

hca = h(isub); isub = isub + 1;
iy = 10;
%pcolor(hca,squeeze(X_xyz(:,iy,:)),squeeze(Z_xyz(:,iy,:)),squeeze(DN_xyz(:,iy,:))')
pcolor(hca,squeeze(DN_xyz(:,iy,:))')

hca = h(isub); isub = isub + 1;
plot(hca,z,dn(Z,lz,R,lr),z,dn(z,lz,0,100*lr),'k--')
clear legends
for ir = 1:numel(r)
  legends{ir,1} = sprintf('r/l_r = %g',r(ir)/lr);
end
%legend(hca,legends,'box','off')
irf_legend(hca,legends,[0.05 0.05])
hca.Title.String = sprintf('lr/lz = %g',lr/lz);


%% Compare different lr
units = irf_units;
dn = @(z,lz,r,lr) units.eps0/units.e*(r.^2./lr.^4-2./lr.^2-1./lz.^2+z.^2./lz.^4).*exp(-z.^2/2./lz - r.^2/2./lr);
dn_xyz = @(x,lx,y,ly,z,lz) units.eps0/units.e*(1./lx.^2+1./ly.^2+1./lz.^2-x.^2./lx.^4-y.^2./ly.^4-z.^2./lz.^4).*exp(-x.^2/2./lx-y.^2/2./ly-z.^2/2./lz);

% Set up figure
nrows = 2;
ncols = 3;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

lz = 1;
ly = 1;
lx = 1;
lr_vec = [0.8 1 1.2 2 5 100]*lz;

for irr = 1:numel(lr_vec)
  lr = lr_vec(irr);
  z = linspace(-10*lz,10*lz,101);
  x = z;
  y = z;
  r = [0 0.5 1 1.5 2]*lr;
  [Z,R] = meshgrid(z,r);
  DN = dn(Z,lz,R,lr);
  [X_xyz,Y_xyz,Z_xyz] = meshgrid(x,y,z);
  DN_xyz = dn_xyz(X_xyz,lx,Y_xyz,ly,Z_xyz,lz);

  DN_xz = sum(DN_xyz,2);
  DN_z = sum(DN_xz,1);


  hca = h(isub); isub = isub + 1;
  plot(hca,z,dn(Z,lz,R,lr),z,dn(z,lz,0,100*lr),'k--')
  clear legends
  for ir = 1:numel(r)
    legends{ir,1} = sprintf('r/l_r = %g',r(ir)/lr);
  end
  irf_legend(hca,legends,[0.05 0.05])
  hca.Title.String = sprintf('lr/lz = %g',lr/lz);
end