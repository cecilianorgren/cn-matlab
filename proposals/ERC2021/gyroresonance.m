% Gyroresonance:
% Relativistic. In relativistic limit, v is limited by c, therefore, it
% becomes important to take into account the distributions between vpar and
% vperp.
% omega - kpar*vpar + n*omega_ce/gamma, where n = 0,+-1,+-2,... and gamma = 1/sqrt(1-v^2/c^2)
% Parallel resonance
% omega - kpar*vpar - n*omega_ce, where n = 0,+-1,+-2,... n = 0 is Landau resonance
units = irf_units;
n = 1e-6; % m-3
vt = 1000e3; % m/s
vtpar = vt;
vtperp = 2*vt;
vtx = vtperp;
vty = vtperp;
vtz = vtpar;
vdx = 0; % m/s
vdz = 0; % m/s
vdy = 0; % m/s
nv = 51;
v = 5*vt*linspace(-1,1,nv);

[VX,VY,VZ] = ndgrid(v,v,v);
VTOT = sqrt(VX.^2 + VY.^2 + VZ.^2);
B = 10e-9; % T
L = 1000e3; % m/s
F = fmax(VX,VY,VZ,n,vdx,vdy,vdz,vtx,vty,vtz);
VPERP = sqrt(VX.^2 + VY.^2);
PA = atan2d(VPERP,VZ);
PA(VTOT>3*vt) = NaN;

oce = units.e*B/units.me;
vph = 1000e3; % = omega/kpar, m/s
k = 2*pi/1000; % s/m

% Resonance condition
RES = abs(VZ-vph)-oce/k;

hca = subplot(1,1,1);
pcolor(hca,v*1e-3,v*1e-3,squeeze((RES(:,ceil(nv/2),:)))')
shading(hca,'flat')
hold(hca,'on')
contour(hca,v*1e-3,v*1e-3,squeeze((F(:,ceil(nv/2),:)))','color',[0 0 0])
hold(hca,'off')
hca.XLabel.String = 'v_{\perp,1} (km/s)';
hca.YLabel.String = 'v_{||} (km/s)';