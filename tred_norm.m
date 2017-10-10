function out = tred_norm(in,field,omega_pi)
% TRED_NORM(in,field,omega_pi)
% Takes TRED46 outputs and converts them to typical space physichs unities,
% i.e. nT, mV/m, km/s, cc
% 
% in - tred46 quantity
% field - define quantity; 'B','E','Ji','Je','Ni','Ne','v','l','t'
% omega_pi - ion plasma frequency of the environement you want to compare
%            the simulation to

% Define constants
c=2.99792458e10; % cm/s
e=4.80320427e-10; % Fr
mp=1.6726e-24; % g

% Define normalizations
switch field
    case 'B' % B
        prim=omega_pi*c*mp/e;
        out=in*prim; % cgs, Gauss
        out=out*1e-4; % SI, T
        out=out*1e9; % nT
    case 'E' % E
        prim=omega_pi*c*mp/e;
        out=in*prim; % cgs, statV/cm
        out=out*1e-6*c; % SI, V/m
        out=out*1e3; % mV/m
    case 'v' % v
        prim=c;
        out=in*prim; % cgs, cm/s
        out=out*1e-2; % SI, m/s
        out=out*1e-3; % km/s
    case 't' % t
        % 1 dt=0.125 omepa_pi^-1
        prim=0.125/omega_pi;
        out=in*prim; % cgs = SI = s
    case 'l' % l
        prim=c/omega_pi;
        out=in*prim; % cgs, cm
        out=out*1e-2; % SI, m
        out=out*1e-3; % km
    case 'Ji' % to v
        % Current density. j=nev ? 
        % 
        prim=c;
        out=in*prim; % cgs, cm/s
        out=out*1e-2; % SI, m/s
        out=out*1e-3; % km/s
    case 'Je'
        % Gl?m inte bort normaliseringen f?r densiteten.
        prim=c;
        out=in*prim; % cgs, cm/s
        out=out*1e-2; % SI, m/s
        out=out*1e-3; % km/s       
    case 'Ni' % seems high
        % TRED46 delivers charge density, mult. by 4*pi and sign species 
        % charge to get density.
        out=in*4*pi*1; % cgs, cm^-3
        
    case 'Ne'
        % TRED46 delivers charge density        
        out=in*4*pi*(-1); % cgs, cm^-6
    otherwise
        out=in;
end
