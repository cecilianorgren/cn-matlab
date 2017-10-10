function  x_res = eom(t,x_vect)
G = 6.67384e-11; % m^3 kg^-1 s^-2 (N m^-2 kg^-2)
ME = 5.9722e24; % kg

x = x_vect(1);
y = x_vect(2);
vx = x_vect(3);
vy = x_vect(4);

x_res = zeros(4,1);

x_res(1) = vx; % dx/dt = vx;
x_res(2) = vy; % dy/dt = vy;
x_res(3) = -sign(x)*G*ME/(x^2+y^2);
x_res(4) = -sign(y)*G*ME/(x^2+y^2);
                