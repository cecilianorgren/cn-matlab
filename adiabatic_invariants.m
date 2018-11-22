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