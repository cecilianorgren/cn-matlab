function  x_res = eom(t,rv,m,q,Ex,Ey,Ez,Bx,By,Bz)

x = rv(1);
y = rv(2);
z = rv(3);
vx = rv(4);
vy = rv(5);
vz = rv(6);

Ex = Ex(x,y,z);
Ey = Ey(x,y,z);
Ez = Ez(x,y,z);
Bx = Bx(x,y,z);
By = By(x,y,z);
Bz = Bz(x,y,z);

% Equations to be solved
x_res = zeros(6,1);
x_res(1) = vx; % dx/dt = vx;
x_res(2) = vy; % dy/dt = vy;
x_res(3) = vz; % dz/dt = vz;
x_res(4) = (q/m)*(Ex + vy*Bz - vz*By);
x_res(5) = (q/m)*(Ey + vz*Bx - vx*Bz);
x_res(6) = (q/m)*(Ez + vx*By - vy*Bx);                                              

