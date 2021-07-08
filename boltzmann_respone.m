% boltzmann_respone
n0 = 1;
T = 1;
e = 1;
x = linspace(-3,3,100);
phi = @(x) exp(-x.^2);
n = @(phi,T,q) n0*exp(-q*phi/T);

hca = subplot(2,1,1);
plot(hca,x,phi(x))

hca = subplot(2,1,2);
plot(hca,x,n(phi(x),T,-e),x,n(phi(x),T,+e))
legend('q = -1','q = +1')