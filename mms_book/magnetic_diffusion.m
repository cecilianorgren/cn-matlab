B0 = 1;
L = 1;
fBx = @(z) B0*tanh(z/L);

zv = 5*L*linspace(-1,1,500);
dz = z(2) - z(1);
%zv_d1 = zvdiff(zv)

Bx = fBx(zv);
d2Bx = del2(Bx,dz);

h = setup_subplots(2,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,zv,Bx)
hca.XLabel.String = 'z/L';
hca.YLabel.String = 'B_x';

hca = h(isub); isub = isub + 1;
plot(hca,zv,d2Bx)
hca.XLabel.String = 'z/L';
hca.YLabel.String = '-\nabla^2 B_x';