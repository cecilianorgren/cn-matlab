

B0 = 1;
L = 1;

Bx_Harris = @(z) -B0*tanh(z/L);
Jy_Harris = @(z) -B0*(1 - tanh(z/L).^2)/L;


z = linspace(-L*4,L*4,100);
  
hca = subplot(1,1,1);
plot(hca,z,Bx_Harris(z),z,Jy_Harris(z))
hca.XLabel.String = 'z/L';
hca.XGrid = 'on',
hca.YGrid = 'on',
%%
B0 = 1e-9;
L = 1;

Bx_Harris = @(z) -B0*tanh(z/L);
%Bx_HarrisL = @(z,L) -B0*tanh(z/L);
Jy_Harris = @(z) -B0*1e-9*(1 - tanh(z/L).^2)/units.mu0/(L*1e3)*1e9;


z = linspace(-L*4,L*4,100);

h
plot(z,Bx_Harris(z),z,Jy_Harris(z))


