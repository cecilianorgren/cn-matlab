% Compare expressions by chen2005 to data.
% How valid is this approach? They dont take into account a finite drift of
% the EH. And exactly where the EH is in the electron distribution plays a
% large role of how many particles are available to trap. For example, if
% the EH is at velocities far out on the distribution, less electron are
% availbale to trap. Should less electrons indicate less possible trapping
% and smaller allowed potential? Why is the potential amplitude unrelated
% to the plasma density? Maybe change in density adjusts the length scales,
% which in turn are normalized?

% Shouldn't only the phase space density at the EH separatrix be relevant?

% For data, run figure: mms_20170706_005403.figures_selfcontained

% Check basic expressions of Chen.
Te = 1;
Ti = 1*Te; % put to low value to exclude ion dynamics
q = 1;

t = q*Te/Ti;
psi = phi/Te;

fun_Fp = @(psi,t) pi.*exp(psi).*(1-erf(sqrt(psi))) - exp(-t.*psi).*sqrt(t).*erfi(sqrt(t.*psi))
