% Start from vector and scalar potential, define Lagrangian, Hamiltonian,
% and canonical momentum. use syms to derive equations of motions. 

% Note, E = -grad(phi) - dA/dt (partial). We must choose from which
% potential the electric field is derived, i.e. phi or A. It cannot be
% both (or at least not fully both).

syms L H K V x y z vx vy vz t px_of_v_A py pz Ax Ay Az phi a b m q Erec
syms K_p(px,py,pz)
syms K_v(vx,vy,vz)

% Constants
units = irf_units;
%m = units.me;
%q = -units.e;

% Current sheet structure
%a = 5;
%b = 1;

% Reconnection electric field
%Erec = 1;

% Define potentials. They are gauge dependent. Chose gauge below.
gauge = 'weyl';
switch gauge
  case 'weyl'
    Ax_fun = 0;
    Ay(x,y,z,t) = x^2/a^2 - z^2/b^2 - Erec*t; % Ay'  = Ay + int_0^2(Ey)dt, here we assume a constant Erec, so it just becomes Erec*t.
    Az_fun = 0;
    phi_fun = 0;
  case 'potential'
    Ax_fun = 0;
    Ay_fun = x^2/a^2 - z^2/b^2;
    Az_fun = 0;
    phi_fun = -y*Erec;
end

%Ay = x^2/a^2 - y^2/b^2 - Erec*t;

% Canonical momenta
px_fun = px == m*vx + q*Ax;
py_fun = py == m*vy + q*Ay;
pz_fun = pz == m*vz + q*Az;

%px_fun = m*vx + q*Ax_fun;
%py_fun = m*vy + q*Ay_fun;
%pz_fun = m*vz + q*Az_fun;

vx_fun = isolate(px_fun,vx);
vy_fun = isolate(py_fun,vy);
vz_fun = isolate(pz_fun,vz);

vx_of_p_A = solve(px,vx);
vy_of_p_A = solve(py,vy);
vz_of_p_A = solve(pz,vz);

% Kinetic energy
K = (m/2)*(vx_fun^2 + vy_fun^2 + vz_fun^2); % K as excplicit as possible, i.e. in terms of x, y, x, a, b, m, q
K_p = (m/2)*(vx_of_p_A^2 + vy_of_p_A^2 + vz_of_p_A^2);
K_v = (m/2)*(vx^2 + vy^2 + vz^2);

% Potential energy
V = q*phi;

% Lagrangian
L = K - V;
L_p = K_p - V;
L_v = K_v - V;

% Hamiltonian
H = K + V;
H_p = K_p + V;
H_v = K_v + V;

% Equations of motions

%vx_ = gradient(H_p,px)