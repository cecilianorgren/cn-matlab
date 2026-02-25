% density_depletion
units=irf_units;
lz=5; % km
lr=12; % km
phi0=200; % V
phi = @(r,z) phi0*exp(-0.5*r.^2/lr/lr-0.5*z.^2/lz/lz);
ld=2.1;
% solve poissons equation del2phi = rho/eps0 to get rho0
% cylindrical coordinates


% cartesian, 1D hole
r=0;
% cc
% but this is mass density not number density, right?
% rho = @(z,r,lz,lr) -units.eps0*phi0*((z.^2/lz/lz-1)/lz/lz+(r.^2/lr/lr-2)/lr/lr).*exp(-0.5*z.^2/lz/lz-0.5*r.^2/lr/lr);
% so to get relative density changes, we divide by the background mass
% density rho_bg = units.me*0.04*1e-6;
nbg_cc = 0.04;
rho_bg = nbg_cc*1e6*units.me;
rho_bg=1/nbg_cc/1e6/10000;
rho = @(z,r,lz,lr) -units.eps0*phi0*((z.^2/lz/lz-1)/lz/lz+(r.^2/lr/lr-2)/lr/lr).*exp(-0.5*z.^2/lz/lz-0.5*r.^2/lr/lr);
z=linspace(-4*lz,4*lz,100);
plot(z/ld,rho(z,0,lz,lz)/rho_bg,...
     z/ld,rho(z,0,lz,lr)/rho_bg,...
     z/ld,rho(z,0,lz,lr*10)/rho_bg)
legend('lr=lz','lr=2.4lz','lr=10lz')