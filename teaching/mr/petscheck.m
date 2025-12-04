B1x = @(x,z) log((4/pi).*((sin(x).^2+sinh(z).^2)./(cos(x).^2+sinh(z).^2)).^0.5);
B1z = @(x,z) atan(sin(2*x)./sinh(2*z));


xx = linspace(-1,1,100);
plot(xx,B1x(xx,0),xx,B1z(xx,0))

