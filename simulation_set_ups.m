% simulations set ups
% Hesse 2016b
alpha = 14.87;
l = 0.5;
Bx = @(z) (0.5 + tanh(z/l))*cosd(alpha) + sind(alpha);
By = @(z) -(0.5 + tanh(z/l))*sind(alpha) + sind(alpha);

z = linspace(-3,3,100);

plot(z,Bx(z),z,By(z))
legend('B_x','B_y')
title('Magnetic field, Hesse et al. 2016b')
xlabel('x')
ylabel('B')

