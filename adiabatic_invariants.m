% adiabatic_invariants

% Isotropic distribution
np = 1000;
angle = 0  + (180-0)*rand(np,1);
nv = 100;
vmax = 1;
v = 0  + (vmax-0)*rand(np,1);
E = v.^2/2;
phi = 0.03;

vperp0 = v.*sind(angle);
vpar0 = v.*cosd(angle);
B0 = 1;

% mu*B = Eperp;
mu = vperp0.^2/2/B0;

vperp = @(B) sqrt(2*mu*B);
vpar = @(B) sign(vpar0).*sqrt(2*E-2*mu*B+2*phi);
vpar_trapped = @(B) sign(vpar0).*sqrt(2*E-2*mu*B+2*phi);
vpar_passing = @(B) sign(vpar0).*sqrt(2*E-2*mu*B+2*phi);


B1 = 0.5*B0;

figure(200)
nrows = 2;
ncols = 2;
npanels = ncols*nrows;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

hca = h(isub); isub = isub + 1;
scatter(hca,vpar0,vperp0)
axis(hca,'equal')
hca.YLim = [0 vmax];
hca.XLim = [-vmax vmax];

hca = h(isub); isub = isub + 1;
scatter(hca,vpar(B1),vperp(B1))
axis(hca,'equal')
hca.YLim = [0 vmax];
hca.XLim = [-vmax vmax];

hca = h(isub); isub = isub + 1;
scatter(hca,vpar0,vpar(B1))
axis(hca,'equal')
hca.YLim = [0 vmax];
hca.XLim = [-vmax vmax];

hca = h(isub); isub = isub + 1;
scatter(hca,vperp0,vperp(B1))
axis(hca,'equal')
%hca.YLim = [0 vmax];
%hca.XLim = [-vmax vmax];


%% Egedal 2008, p. 9-10-
% distribution function
ninfty = 1;
vtpar = 1;
vtperp = 0.5;
f = @(vpar,vperp1,vperp2) ninfty/sqrt(2*pi^3*vtpar*vtperp^2)*exp(-0.5*(vpar.^2/vtpar^2+vperp1.^2/vtperp^2+vperp2.^2/vtperp^2));

nvperp = 110; 
nvpar = 110;
vmax = 3*vtpar;
vperp = 3*linspace(-vmax,vmax,nvperp);
vperp1 = vperp;
vperp2 = vperp;
vpar = 2*linspace(-vmax,vmax,nvpar);
dvperp = vperp1(2)-vperp1(1);
dvpar = vpar(2)-vpar(1);
d3v = dvpar*dvperp^2;

[VPAR,VPERP1,VPERP2] = meshgrid(vpar,vperp1,vperp2);
F = f(VPAR,VPERP1,VPERP2);
nint = sum(F(:))*d3v;

hca = subplot(1,1,1);
pcolor(hca,squeeze(VPAR(:,:,1)),squeeze(VPERP1(:,:,1)),sum(F,3));
hca.Title.String = sprintf('ninfty = %g, nint = %g',ninfty,nint);
colorbar

