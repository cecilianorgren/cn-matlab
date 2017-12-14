%% Set up
units = irf_units;

% Plasma properties
n = 0.06e6;
ne = n;
ni = n;
Te = 200:200:3000; Te_K = Te*units.eV/units.kB;  % use perpendicular temperature -> diamagnetic drift
Ti = 7000; Ti_K = Ti*units.eV/units.kB; % use perpendicular temperature -> diamagnetic drift
B = 10e-9; 

% Physical constants
qe = 1.602e-19;
me = 9.109e-31;
mi = 1.673e-27;
mime = 1836;
eps = 8.854e-12;
kb = 1.381e-23;

% Physical parameters
vti = sqrt(2*qe*Ti/mi);
vte = sqrt(2*qe*Te/me);
wpe = sqrt(ne*qe^2/(me*eps)); % Hz
wpi = sqrt(ni*qe^2/(mi*eps)); % Hz
wce = qe*B/me;
wci = qe*B/mi;
roe = vte/wce;
roi = vti/wci;
Ld = sqrt(vte.^2)/wpe/sqrt(2)*1e-3;
wlh = sqrt(wce*wci);
PB = B^2/2/units.mu0;
Pi = Ti_K*units.kB*n;
Pe = Te_K*units.kB*n;
betai = Pi/PB;
betae = Pe/PB;
beta = betai + betae;


% Gradient length scales, m
% E = -50*1e-3; % E (together with Ln/LT) is used to get Ln and other parameters
% LnLT = 0.1; % Ln/LT, Ln generally smaller then LT
vE = -E/B;
vdi = -vE;
vde = -Te/Ti*vdi;

Ln = -(1+LnLT)*(Te_K*units.kB/units.e/B)./vde;
LT = Ln/LnLT;
LB = -2./(beta./Ln + betae./LT); % Davidson1975 eq. 13: (1/LB) = 0.5*beta*(1/Ln) - 0.5*betae*(1/LT)

vB = -1./LB.*vte.^2/2/wce;
vn = -1./Ln.*vte.^2/2/wce;
vT = -1./LT.*vte.^2/2/wce;    
plot(Te,[vn' vT' vB'],...
     Te,vB-vn,Te,vph_allmax,...
     Te,vph_allmax-vE,...
     Te,repmat(vti,1,numel(Te)),...
     Te,(vB-vn)./(vph_allmax-vE) ...
   ); 
legend('vn','vT','vB','vB-vn','vph','vph-vE','vti','(vB-vn)/(vph-vE)')
