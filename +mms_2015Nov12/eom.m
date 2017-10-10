function  x_res = eom(t,x_vect,a,B0,d,E0,eps)
% physical constants
e = 1.6022e-19;
me = 9.10939999e-31;
% field parameters

x = x_vect(1);
y = x_vect(2);
z = x_vect(3);
vx = x_vect(4);
vy = x_vect(5);
vz = x_vect(6);

Bx = -B0*tanh(z*pi/d).*(1-exp(-z.^2*5/(d^2)))-0*B0/6;
By = 1*0.5*B0-3*0.5*B0*sin(4/3*pi/d*z).*exp(-z.^2*2/(d^2)).*(1-exp(-z.^2*5/(d^2)));%0.5*B0-2*0.5*B0*sin(1*pi/d*z).*exp(-z.^2*5/(d^2));%0.5*B0-0.5*B0*sin(pi/d*z);     
Bz = B0*eps;
Ex = 0;
Ey = 0;
Ez = E0*sin(pi/d*z)*1;

x_res = zeros(6,1);
x_res(1) = vx; % dx/dt = vx;
x_res(2) = vy; % dy/dt = vy;
x_res(3) = vz; % dz/dt = vz;
x_res(4) = (-e/me)*(Ex + vy*Bz - vz*By);
x_res(5) = (-e/me)*(Ey + vz*Bx - vx*Bz);
x_res(6) = (-e/me)*(Ez + vx*By - vy*Bx);                                              
                