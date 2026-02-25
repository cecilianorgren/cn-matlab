B0 = 1;
L = 1;
BG = 0.5;
Bz = @(z) B0*tanh(z/L);
By = @(z) BG*B0*z./z;

z = linspace(-2*L,2*L,100);

subplot(2,1,1)
plot(z,[Bz(z);By(z)])

subplot(2,1,2)
plot(z,[Bz(z);By(z)])

