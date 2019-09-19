% Comparison to Mozer 2018 model

% They show v inside as a function of v outside, so basically v(v0).

units = irf_units;

fun_v = @(v0,vph,phi) vph + sign(v0-vph).*((v0-vph).^2 + 2*units.e*phi/units.me).^0.5;

%V0(iabove) = VPH(iabove) + ((V(iabove)-VPH(iabove)).^2 + 2*q*(PHI(iabove))/m).^0.5; % with proper substititions, same as Schamel free streaming
%V0(ibelow) = VPH(ibelow) - ((V(ibelow)-VPH(ibelow)).^2 + 2*q*(PHI(ibelow))/m).^0.5;
%V0(itrap) = NaN;

vph = -10000*1e3; % m/s
phi = 500; % V

v_max = 3e7; % m/s
v_vec = linspace(-v_max,v_max,1000); % m/s

hca = subplot(1,1,1);

plot(hca,v_vec,fun_v(v_vec,vph,phi))
hca.YLim = v_max*[-1 1];
hca.XGrid = 'on';
hca.YGrid = 'on';

