function [t x y] = orbit(r_per,r_ap,dt,Ttot)
G = 6.67384e-11; % m^3 kg^-1 s^-2 (N m^-2 kg^-2)
ME = 5.9722e24; % kg
RE = 6371*1e3; % m

%r_per = (3+1)*RE; % m, perigee
%r_ap = (15+1)*RE; % m, apogee
a = (r_ap + r_per)/2; % m, semi-major axis

mu = G*ME; % m^3 s^-2
e = (r_ap - r_per)/(r_ap + r_per); % eccentricity ?
V_ap = sqrt( (1-e)*mu/(1+e)/a ); % m/s, velocity at apogee
V_peri = sqrt( (1+e)*mu/(1-e)/a ); % m/s, velocity at perigee

n = sqrt(mu/a^3);
T = 2*pi/n; % orbital period
%nt = 1000;
%t = linspace(0,T,nt);
%dt = 120; % s, 120s=2min
t = 0:dt:T;
nt=numel(t)
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
tstop = Ttot;60*60*24*365;
nT = ceil(tstop/T);
rr = repmat(r,1,nT);
ff = repmat(f,1,nT);
tt = linspace(0,nT*T,nT*nt);
ffplus = tt*2*pi/(60*60*24*365);
x = rr.*cos(ff+ffplus);
y = rr.*sin(ff+ffplus);
