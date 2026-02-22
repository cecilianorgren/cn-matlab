hz = @(Pz,x,z) 0.5*(Pz.^2 + (x - z.^2/2).^2);

z = linspace(-3,3,100);
pz = linspace(-3,3,100);
[Z,PZ] = ndgrid(z,pz);

HZ = hz(PZ,02,Z);

contour(Z,PZ,HZ)

%%
ip = 2;
tt = particles{ip}.t;
zz = particles{ip}.z;
yy = particles{ip}.y;
vvz = particles{ip}.vz;

nrows = 3;
ncols = 1;
isub = 1;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,zz,yy) 
hca.XLabel.String = 'z';
hca.YLabel.String = 'y';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,zz,vvz) 
hca.XLabel.String = 'z';
hca.YLabel.String = 'v_z';

hca = subplot(nrows,ncols,isubresults); isub = isub + 1;
Iz = cumtrapz(zz,vvz);
plot(hca,tt,Iz) 
hca.XLabel.String = 't';
hca.YLabel.String = 'I_z';

%hca = subplot(nrows,ncols,isub); isub = isub + 1;
%Iz = cumtrapz(zz,vvz);
%plot(hca,Iz) 