data_tmp = load('/Users/cno062/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat');
obs_eh_properties = data_tmp.EH_properties;
obs_lpp = obs_eh_properties.Lpp; % peak to peak length
obs_potential = obs_eh_properties.potential;
obs_potential_max = obs_eh_properties.max_potential;
obs_velocity = obs_eh_properties.vel;
obs_neh = numel(velocity);


phi = @(x,l,phi0) phi0*exp(-x.^2/2./l.^2);
neni = @(x,l,phi0) -units.eps0/units.e*(1/l^2)(1 + (x/l).^2)*phi0.*exp(-x.^2/2/l^2);

ieh = 1;

x = linspace(-3*obs_lpp(ieh),3*obs_lpp(ieh),100);

hca = subplot(1,1,1);
plot(hca,x,phi(x,obs_lpp(ieh),obs_lpp(ieh)));

% phi = phi0*exp(-x^2/2/l^2), 
% dphi/dx = (-2*x/2/l^2)*phi0*exp(-x^2/2/l^2)
% d2phi/dx2 = (-1/l^2)*phi0*exp(-x^2/2/l^2) + (-2*x/2/l^2)^2*phi0*exp(-x^2/2/l^2) 
%           = [(-1/l^2) + (-x/l^2)^2]*phi0*exp(-x^2/2/l^2) 
%           = -(1/l^2)[1 + (x/l)^2]*phi0*exp(-x^2/2/l^2) 
%           == -(e/eps0)*(ni-ne)
% @ x = 0: 
% d2phi/dx2 = -(1/l^2)*phi0 = -(e/eps0)*(ni-ne)
%   => (ni-ne) = eps0/e*(1/l^2)*phi0
% potential from charge density
dn = units.eps0/units.e*obs_potential_max./(obs_lpp*1e3)*1e-6; % cc