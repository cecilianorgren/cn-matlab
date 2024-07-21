%% Trying to get some function that can explain ion distribution in (vL,vM)

% H = (m/2)*(vx^2 + vy^2 + vx^2) + q*phi
%   = (m/2)*[(px-q*Ax)^2 + (py-q*Ay)^2 + (pz-q*Az)^2]  + q*phi
% Ax = -z^2/2*By0/Lz
% phi = -Ez0*z/Lz

Lz = 1;
By0 = 1;
Ez0 = 3;
q = 1;
m = 1;


Ax = @(z) z.^2/2*By0/Lz;
phi = @(z) Ez0*z.^2/2/Lz;

By = @(z) z*By0/Lz; % By = (d/dz)Ax;
Ez = @(z) -Ez0*z/Lz;

z0 = 2;
vx0 = 0;

vx = @(z) vx0 + (By0/Lz)*(z.^2-z0.^2);

dUpot = @(z) 2*(phi(z)-phi(z0));

vz2 = @(z) -dUpot(z) - vx(z).^2;

vz2(0)

zvec = 2*linspace(-Lz,Lz,100);

nrows = 4;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,zvec,Ax(zvec),zvec,By(zvec))
hca.XLabel.String = 'z (...)';
irf_legend(hca,{'A_x','B_y'},[0.02 0.98])


hca = h(isub); isub = isub + 1;
plot(hca,zvec,phi(zvec),zvec,Ez(zvec))
hca.XLabel.String = 'z (...)';
irf_legend(hca,{'\phi','E_z'},[0.02 0.98])


hca = h(isub); isub = isub + 1;
plot(hca,zvec,vx(zvec))
hca.XLabel.String = 'z (...)';
irf_legend(hca,{'vx'},[0.02 0.98])

hca = h(isub); isub = isub + 1;
plot(hca,zvec,vx(zvec).^2,zvec,vz2(zvec),zvec,dUpot(zvec))
hca.XLabel.String = 'z (...)';
irf_legend(hca,{'vx^2','vz^2','dU_{pot}'},[0.02 0.98])
