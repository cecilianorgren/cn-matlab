%% Trying to get some function that can explain ion distribution in (vL,vM)

% H = (m/2)*(vx^2 + vy^2 + vx^2) + q*phi
%   = (m/2)*[(px-q*Ax)^2 + (py-q*Ay)^2 + (pz-q*Az)^2]  + q*phi
% Ax = -z^2/2*By0/Lz
% phi = -Ez0*z/Lz

Lz = 1;
By0 = 1;
Ez0 = 1;

Ax = @(z) -z.^2/2*By0/Lz;
phi = @(z) Ez0*z.^2/2/Lz;
Ez = @(z) -Ez0*z/Lz;

By = @(z) z*By0/Lz; % By = -dAx/dz;


zvec = linspace(-3*Lz,3*Lz,100);

nrows = 3;
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


