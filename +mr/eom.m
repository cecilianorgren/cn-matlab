function  x_res = eom(t,x_vect,a,b,d,eps)
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

Bx = b*z/d;
By = 0;
Bz = -b*eps+0;
%Bz = b*x/d/15;
Ex = 0;
Ey = a;
Ez = 0;

x_res = zeros(6,1);
x_res(1) = vx; % dx/dt = vx;
x_res(2) = vy; % dy/dt = vy;
x_res(3) = vz; % dz/dt = vz;
x_res(4) = (-e/me)*(Ex + vy*Bz - vz*By);
x_res(5) = (-e/me)*(Ey + vz*Bx - vx*Bz);
x_res(6) = (-e/me)*(Ez + vx*By - vy*Bx);                                              
                