%
units = irf_units;
U = @(v,vphi,phi) 0.5*units.me.*(v-vphi).^2-units.e*phi;
integrand = @(U,a) exp(-U)./sqrt(U-a);
indefinite_integral = @(U,a) sqrt(pi).*exp(a).*(erf(sqrt(U+a))-1);
definite_integral = @(U1,U2,a) indefinite_integral(U2,a) - indefinite_integral(U1,a);

x1 = 0;
x2 = 100;
nx = 200;
x = linspace(x1,x2,nx);
phi_max = 1000/220;%1800/220; 
phi_L = (x2-x1)/20;
phi_x = (x2-x1)/2;
phi = @(x) phi_max*(1+tanh((x-phi_x)/phi_L));


h = setup_subplots(3,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,x,phi(x))
hca.XLabel.String = 'x';
hca.YLabel.String = 'phi(x)';

hca = h(isub); isub = isub + 1;
plot(hca,x,definite_integral(0,Inf,phi(x))./definite_integral(0,Inf,phi(x(1))))
hca.XLabel.String = 'x';
hca.YLabel.String = 'n(phi(x))';

hca = h(isub); isub = isub + 1;
plot(hca,phi(x),definite_integral(0,Inf,phi(x))./definite_integral(0,Inf,phi(x(1))))
hca.XLabel.String = 'phi(x)';
hca.YLabel.String = 'n(phi(x))';


