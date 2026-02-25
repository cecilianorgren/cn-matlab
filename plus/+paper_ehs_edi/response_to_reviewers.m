% Figures for respone to reviewers


%% Non-linear relation between dn and phi
f_phi = @(x,lx,phi0) phi0.*exp(-x.^2./(lx.^2));
f_E = @(x,lx,phi0) phi0.*exp(-x.^2./(lx.^2));

lx = 1;
x = linspace(-3*lx,3*lx,100);

plot(x,f_phi(x,lx,1))

%% Non-linear relation between dn and phi, symbolic expressions
units = irf_units;
num_lx = 1;
%phimax = 1;
syms x v phimax lx %phi(x) f0(v) ff(x,v) v0(x,v)
f_phi(x,phimax,lx) = phimax*exp(-x.^2/2/lx.^2);
f_E = -diff(f_phi,x);
f_dn = -diff(units.eps0/units.e*f_phi,x,2); % dn = -diff(E,x); % the same, dn = ne-ni

vecx = linspace(-3,3,100);
vecphimax = 0:100;
veclx = linspace(0.5,5,100);

h = setup_subplots(4,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,vecx,double(f_phi(vecx,1,1)))
hca.XLabel.String = 'x';
hca.YLabel.String = '\phi';

hca = h(isub); isub = isub + 1;
plot(hca,vecx,double(f_E(vecx,1,1)))
hca.XLabel.String = 'x';
hca.YLabel.String = 'E';

hca = h(isub); isub = isub + 1;
plot(hca,vecx,double(f_dn(vecx,1,1)))
hca.XLabel.String = 'x';
hca.YLabel.String = '\delta n';

hca = h(isub); isub = isub + 1;
scatter(hca,veclx,double(f_dn(0,1,veclx)))
hca.XLabel.String = 'l_x';
hca.YLabel.String = '\delta n';
%%
% get from paper_fig1.m
% phi1,phi2,phi3,phi4
dn1 = dn_E1par/v_for_density_scaling*1e-6/units_scaling; dn1 = dn1.resample(phi1);
dn2 = dn_E2par/v_for_density_scaling*1e-6/units_scaling; dn2 = dn2.resample(phi2);
dn3 = dn_E3par/v_for_density_scaling*1e-6/units_scaling; dn3 = dn3.resample(phi3);
dn4 = dn_E4par/v_for_density_scaling*1e-6/units_scaling; dn4 = dn4.resample(phi4);

hca = subplot(1,1,1);
scatter(hca,phi1.data,dn1.data)
hold(hca,'on')
scatter(hca,phi2.data,dn2.data)
scatter(hca,phi3.data,dn3.data)
scatter(hca,phi4.data,dn4.data)
hold(hca,'off')