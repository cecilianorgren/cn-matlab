
lambda = 1;
k = 2*pi/lambda;

nx = 100;
ny = 100;
x = linspace(0,3*lambda,nx);
y = linspace(0,3*lambda,ny);

[X,Y] = ndgrid(x,y);

B0_ = 1;
B0 = @(X,Y,B0_) X*0+B0_;

angle = 45;
kx = k*sind(angle);
ky = k*cosd(angle);

B1_ = 0.1*B0;
B1 = B1_;