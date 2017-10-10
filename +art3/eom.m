function  x_res = EquationOfMotion(t,x_vect,E)
% linear acceleration from constant electric field
e = 1.6022e-19;
me = 9.10939999e-31;

x = x_vect(1);
vx = x_vect(2);

x_res = zeros(2,1);
x_res(1) = vx; % dx/dt = vx;
x_res(2) = (-e/me)*E; % dvx/dt = ax
                