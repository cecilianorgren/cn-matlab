% isosurface

[y,x,z] = ndgrid(linspace(-.75,.75,100),linspace(-.1,2.1,100),linspace(-.2,.2,100));
f = (x.*(x-1).^2.*(x-2) + y.^2).^2 + z.^2;
cla
isosurface(x,y,z,f,.01);
view(3);
camlight
axis equal

%% isosurface
syms x y z f(x,y,z)
r = [x,y,z];
lx = 1;
ly = 2;
lz = 3;
f0 = 10;

[Y,X,Z] = ndgrid(5*linspace(-lx,lx,9),5*linspace(-ly,ly,10),5*linspace(-lz,lz,11));

f(r) = f0*exp(-0.5*(x/lx).^2 -0.5*(y/ly).^2 -0.5*(z/lz).^2);
%symf = symfun(symf,[x,y,z]);
%f = matlabFunction(symf,[x,y,z]);
F = f(X,Y,Z);
cla
isosurface(X,Y,Z,F,1);
view(3);
%camlight
axis equal
xlabel('x')
ylabel('y')
zlabel('z')