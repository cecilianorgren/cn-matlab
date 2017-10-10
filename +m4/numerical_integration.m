G = 6.67384e-11; % m^3 kg^-1 s^-2 (N m^-2 kg^-2)
ME = 5.9722e24; % kg
RE = 6371*1e3; % m

r_per = (3+1)*RE; % m, perigee
r_ap = (15+1)*RE; % m, apogee
a = (r_ap + r_per)/2; % m, semi-major axis

mu = G*ME; % m^3 s^-2
e = (r_ap - r_per)/(r_ap + r_per); % eccentricity ?
V_ap = sqrt( (1-e)*mu/(1+e)/a ); % m/s, velocity at apogee
V_peri = sqrt( (1+e)*mu/(1-e)/a ); % m/s, velocity at perigee

n = sqrt(mu/a^3);
T = 2*pi/n; % orbital period

% integrate numerically, NOPE, it's to sensitive
% some sun inertial coordinate system, x initially towards sun.
x0 = r_per; % m
y0 = 0;
z0 = 0;
vx0 = 0; % m/s
vy0 = V_peri+5.8;
vz0 = 0;
x_init = [x0;y0;vx0;vy0]; % m, m/s
T = 24*60*60; % one year in seconds

% Integrate trajectory
EoM = @(ttt,xxx) m4.eom(ttt,xxx);
[t,x_sol] = ode45(EoM,[0 T],x_init);        

x = x_sol(:,1);
y = x_sol(:,2);
vx = x_sol(:,3);
vy = x_sol(:,4);

plot(x/RE,y/RE)
