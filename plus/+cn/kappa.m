function out = kappa(v,k,T,n)
% Generates a kappa distribution.
%   f = cn.kappa(v,k,T,n);
%       T - eV
%       n - cc
%       v - km/s
%       k - kappa index

e = 1.6022e-19;
kB = 1.38e-23;
m = 9.1094e-31;

v = v*1e3; % m/s
n = n*1e6; % m^-3
T = T*e/kB; % K

distr=2;
switch distr
    case 1
        w = @(k,T) sqrt((2*k-3)*kB*T/k/m);
        %w = @(k,T) sqrt((2)*kB*T/m);
        f = @(v,k,T,n) n/(pi^(3/2)*(k*w(k,T)^2)^(3/2))...
                        *gamma(k+1)/gamma(k-1/2)...
                        *(1+v.^2/(k)/w(k,T)/w(k,T)).^(-k-1);
    case 2
        w = @(k,T) sqrt((2*k-3)*kB*T/k/m);
        %w = @(k,T) sqrt((2)*kB*T/m);
        f = @(v,k,T,n) n/(2*pi*(k*w(k,T)^2)^(3/2))...
                        *gamma(k+1)/gamma(k-1/2)/gamma(3/2)...
                        *(1+v.^2/(k)/w(k,T)/w(k,T)).^(-k-1);
    case 3
        w = @(k,T) sqrt((2)*kB*T/m);

        w = @(k,T) sqrt((2*k-3)*kB*T/k/m);
        f = @(v,k,T,n) n/(2*pi*w(k,T)^3*k^(3/2))...
                        *gamma(k)/gamma(k-3/2)...
                        *(1+v.^2/(k)/w(k,T)/w(k,T)).^(-k);
end

out = f(v,k,T,n); 