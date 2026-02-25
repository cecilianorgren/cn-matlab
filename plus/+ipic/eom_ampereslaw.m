function  x_res = eom_ampereslaw(t,x_vect)
% linear acceleration from constant electric field
mu0 = 1.2566e-6;
eps0 = 8.8542e-12;

x = x_vect(1);
vx = x_vect(2);

x_res = zeros(2,1);
x_res(1) = vx; % dx/dt = vx;
x_res(2) = (-e/me)*E; % dvx/dt = ax
                