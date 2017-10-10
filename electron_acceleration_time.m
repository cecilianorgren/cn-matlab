% electron acceleration
units=irf_units;
E0 = -0.1*1e-3; % 1 mV/m = V/m
fpe = 2200; % Hz
tstop = 1000/fpe; % s
x0 = 0;
vx0 = 0;
x_init = [x0 vx0];

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