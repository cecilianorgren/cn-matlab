% Background: Electrons that are initially accelerated at the separatrices
% should be more efficiently accelerated by Fermi-process at neutral sheet.
% What is the difference in bulk speed if the poulation prior to reflection
% is hotter/colder. That is, how do the post-reflection bulk speed depend
% on the temperature of the pre-reflected population?
% Here I will try to investigate that.

% First the electrons are accelerated at the separatrices, from v0 to v1. 
% This is a function of an acceleration potential, and relative drift of 
% lobe vs. potential.
% Second, this accelrated population is either thermalized or not.
% The thermalization can be modeled simply by multiplication of the
% temperature. Or can it?
% Third, the electrons are reflected at the neutral sheet due tot he
% reconnection electric field. This reflection/acceleration si a Fermi type
% reflection which we model simply by elastic reflection at moving
% boundary: v_after = -v_before + 2*v_DF.

% Initial distribution is drifting or non-drifting Maxwellian, depending on
% where we start. I.e. out in the lobes, or already accelerated population.


units = irf_units;
fun_f_max = @(n,v,vd,vt) n.*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2);
vDF = 1000e3; % m/s

n1a = 0.01*1e6; % m^-3
T1a = 50; % eV
vd1a = 10000e3; % ms-1
vt1a = sqrt(2*units.e*T1a./units.me); % m/s

n1b = 0.01*1e6; % m^-3
T1b = 300; % eV
vd1b = 0000e3; % ms-1
vt1b = sqrt(2*units.e*T1b./units.me); % m/s

vmin = -40000e3; % m/s
vmax = 40000e3; % m/s
nv = 1000;
v1 = linspace(vmin,vmax,nv);

fun_v_acc = @(v0,vph,phi) vph + sign(v0-vph).*((v0-vph).^2 + 2*units.e*phi/units.me).^0.5;
fun_v2 = @(v1,vDF) -v1 + 2*vDF;
%fun_v_refl = 

%fun_f = @(v) exp()


% Plot
h = setup_subplots(2,1); isub = 1;


if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,v1,fun_v2(v1,vDF).^1)
  hca.XLabel.String = 'v_1 (m/s)';
  hca.YLabel.String = 'v_2 (m/s)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

hca = h(isub); isub = isub + 1;
plot(hca,v1,fun_f_max(n1a,v1,vd1a,vt1a),...
         v1,fun_f_max(n1b,v1,vd1b,vt1b))
       
       
hca = h(isub); isub = isub + 1;
plot(hca,v1,fun_f_max(n1a,fun_v2(v1,vDF),vd1a,vt1a),...         
         v1,fun_f_max(n1b,fun_v2(v1,vDF),vd1b,vt1b))
       
h(2).XLim = h(1).XLim;
       
