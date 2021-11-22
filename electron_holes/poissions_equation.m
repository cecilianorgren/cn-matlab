units = irf_units;
q = -units.e;
m = units.me;
%eps = 1;units.eps0;
syms x y z lx ly lz eps0 phi0 phi(x,y,z,lx,ly,lz,phi0) rho(x,y,z,lx,ly,lz,phi0)
inp = [x,y,z,lx,ly,lz,phi0, eps0];
R = [x y z];
phi(inp) = phi0*exp(-0.5*(x/lx)^2-0.5*(y/ly)^2-0.5*(z/lz)^2);

% E = -grad(phi)
E = -gradient(phi,R);

% div(E) = rho/eps0;
rho = eps0*divergence(E,R);
rho0(lx,ly,lz,phi0,eps0) = rho(0,0,0,lx,ly,lz,eps0,phi0);