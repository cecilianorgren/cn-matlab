% Plasma dispersion function

z=-3:0.1:3;
N=16;
plot(z,real(cef(z,N)),z,imag(cef(z,N)));