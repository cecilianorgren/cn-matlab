% model of electron hole

dx = -2:0.2:2;
[x,y,z] = meshgrid(dx,dx,dx);
v = exp(-x.^2 - y.^2 - z.^2);
[px,py,pz] = gradient(v,.2,.2,.2);
%%
quiver3(x,y,z,px,py,pz), hold off