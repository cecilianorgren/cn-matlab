function  x_res = EquationOfMotion(t,x_vect,B0,lr,lz,phi0,L,rlim,zlim,x0,y0)
%global B0 lr lz phi0 L rlim zlim
e = 1.6022e-19;
me = 9.10939999e-31;

x = x_vect(1);
y = x_vect(2);
z = x_vect(3);
vx = x_vect(4);
vy = x_vect(5);
vz = x_vect(6);

x_res = zeros(6,1);
%[B0,lr,lz,phi0,L,rlim,zlim];
x_res(1) = vx; % dx/dt = vx;
x_res(2) = vy; % dy/dt = vy;
x_res(3) = vz; % dz/dt = vz;
x_res(4) = (-e/me)*((x/(lr*1e3)^2).*phi0.*exp(-0.5*((x-x0)/(lr*1e3)).^2-0.5*((y-y0)/(lr*1e3)).^2-0.5*(z/(lz*1e3)).^2) + vy*B0*1e-9);
x_res(5) = (-e/me)*((y/(lr*1e3)^2).*phi0.*exp(-0.5*((x-x0)/(lr*1e3)).^2-0.5*((y-y0)/(lr*1e3)).^2-0.5*(z/(lz*1e3)).^2) - vx*B0*1e-9);
x_res(6) = (-e/me)*((z/(lz*1e3)^2).*phi0.*exp(-0.5*((x-x0)/(lr*1e3)).^2-0.5*((y-y0)/(lr*1e3)).^2-0.5*(z/(lz*1e3)).^2));                                              
                