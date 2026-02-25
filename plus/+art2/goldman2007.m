%% check how goldman 2007 potential looks like compared to a gaussian
lz=1;
phi0=1;
a=1;
gauss = @(z) phi0*exp(-z.^2/2/lz/lz);
gold = @(z,a) phi0*sech(z*(1./a));

z = linspace(-4,4,300);

%plot(z,gauss(z),'--',z,gold(z,0.8),z,gold(z,1),z,gold(z,2),z,gold(z,3))

as=[0.1 0.5 1 2 3];
plot(z,gauss(z),':',z,gold(z',as))
legend('Gaussian','Goldman')
%legends={'Maxwellian',['Goldman: ' num2str()]}

%% table 1 from goldman2007
