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
nt = 200;
t = linspace(0,T,nt);
tau = 0; % time of perihelion passage
M = n*(t-tau);

funEM = @(E,M) E-e*sin(E)-M;
E = fsolve(@(E) funEM(E,M),M);
r = a*(1-e*cos(E));

f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2)); % true anomaly (angle from center of ellipse)
% add procession of orbit as Earth goes around Sun
% 2*pi rad degrees in 60*60*24*365 seconds
fplus = t*2*pi/(60*60*24*365);


% directly from equation
%theta1 = linspace(0,2*pi, nt);
%r1 = a*(1-e^2)./(1+e*cos(theta1));
%
%polar(f+fplus,r,'r');
%polar(theta1,r1,'g')
%legend('r','r1')
%

% duplicate orbit and add the procession
tstop = 60*60*24*6;
nT = fix(tstop/T);
rr = repmat(r,1,nT);
ff = repmat(f,1,nT);
tt = linspace(0,nT*T,nT*nt);
ffplus = tt*2*pi/(60*60*24*365);
x = rr.*cos(ff+ffplus);
y = rr.*sin(ff+ffplus);
plot(x/RE,y/RE,cosd(0:1:360),sind(0:1:360))
axis equal
%%
x = r.*cos(f+fplus);
y = r.*sin(f+fplus);
plot(x/RE,y/RE)
axis equal

disp(['Orbital period: ' num2str(T/60/60,'%.2f') ' days'])
disp(['Velocity at perigee: ' num2str(V_peri*1e-3,'%.2f') ' km/s'])
disp(['Velocity at apogee: ' num2str(V_ap*1e-3,'%.2f') ' km/s'])

%% bin the time spent
xGrid = -r_ap:0.4:r_ap*RE;
[nzz,binz] = histc(z,zGrid*1e3);
[nrr,binr] = histc(r,rGrid*1e3);
        [arz,brz,crz] = unique([binr binz],'rows');
%% integrate numerically, NOPE, it's to sensitive
if 1
% some sun inertial coordinate system, x initially towards sun.
x0 = r_per; % m
y0 = 0;
z0 = 0;
vx0 = 0; % m/s
vy0 = V_peri+5.8;
vz0 = 0;
x_init = [x0;y0;vx0;vy0]; % m, m/s
T = 1*60*60; % one year in seconds

% Integrate trajectory
EoM = @(ttt,xxx) m4.eom(ttt,xxx);
[t,x_sol] = ode45(EoM,[0 P],x_init);        

x = x_sol(:,1);
y = x_sol(:,2);
vx = x_sol(:,3);
vy = x_sol(:,4);

plot(x/RE,y/RE)
end