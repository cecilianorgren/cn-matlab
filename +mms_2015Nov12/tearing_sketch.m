% tearing_sketch
Lx = 1;
Lz = 2;
phi0 = 1;

sf = @(x,z) x.^2*Lz + phi0*cos(2*2*z/Lx/pi);

x = linspace(0,10*Lx,100);
z = linspace(-Lz,Lz,100);

[X,Z] = meshgrid(x,z);

h = subplot(1,1,1);
[c,ch] = contour(X,Z,sf(Z,X),[0:1:0.1:0.9 1:10]); 
clabel(c,ch)
h.Title.String = sprintf('Topology of tearing mode: sf = %g*x^2 + %d*cos(%g*z)',Lz,phi0,3/Lx/2);
h.XLabel.String = 'x';
h.YLabel.String = 'z';
