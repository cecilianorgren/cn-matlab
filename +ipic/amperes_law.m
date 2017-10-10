% Integrate amperes law with LHS = 0 (which is the case in a 1D simulation)
% and see how the current should vary.
% Ampere's law: rot(B) = mu0*J+mu0*eps0*dEdt
% J is current density A/m^2

% I think it makes sense to integrate the current with the given starting
% values. J = J0, E = E0. Or maybe I should try both.

% Initial current
units = irf_units;
n = 0.06*1e6; 
R = 0.5;
S = 0.5;
% It is only the beam that contributes to the initial current
Te2 = 60;
Te1 = 1600;
vte1 = cn_eV2v(1600,'eV');
J = -n*units.e*S*vte1;

%%
dt = 0.01;
ncycles = 20;
TTpe = dt*ncycles/2/pi*sqrt(units.mp/units.me)*sqrt(0.5)
outputcycle = 2;
ToTpe = dt*outputcycle/2/pi*sqrt(units.mp/units.me)*sqrt(0.5)
%%
% Perform integration
fpe = 2200; % Hz
tstop = 1000/fpe; % s
x_init = [J0 E0];

EoM = @(ttt,xxx) art3.eom(ttt,xxx,E0);
[t,x_sol] = ode45(EoM,[0 tstop],x_init);
x = x_sol(:,1); % m
vx = x_sol(:,2); % m/s
Ex = units.me*vx.^2/2;

subplot(2,1,1)
plot(t*fpe,x/1000/1000) % 10^3 km
xlabel('f_{pe}t')
ylabel('x [10^3 km]')
title('Electron acceleration by constant electric field')

subplot(2,1,2)
plot(t*fpe,vx/1000/1000) % 10^3 km/s
xlabel('f_{pe}t')
ylabel('v [10^3 km/s]')
title('Electron acceleration by constant electric field')

if 0
    plot(t*fpe,Ex) % km/s
    xlabel('f_{pe}t')
    ylabel('E_e')
end
