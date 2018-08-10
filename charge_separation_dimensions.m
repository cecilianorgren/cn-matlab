n0 = 1;
lx = 1;
ly = 10;
lz = 1;

n = @(n0,x,y,z,lx,ly,lz) n0*exp(-x.^2/2/lx^2 - y.^2/2/ly^2 -z.^2/2/lz^2);

x = linspace(-3*lx,3*lx,100); dx = x(2)-x(1);
y = linspace(-3*ly,3*ly,100); dy = y(2)-y(1);
z = linspace(-3*lz,3*lz,100); dz = z(2)-z(1);

[X,Y,Z] = meshgrid(x,y,z);
N = n(n0,X,Y,Z,lx,ly,lz);

Ex = diff(N,1,1)/dx;
Ey = diff(N,1,2)/dy;
Ez = diff(N,1,3)/dz;
E = sqrt(Ex.^2 + Ey.^2 + Ez.^2);