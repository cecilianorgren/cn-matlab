% Physical variables:
B = 25;
n = 0.04;
R = 0.1;
n1= n*(1-R);
n2= n*R;
no = 0;
Te1 = 1600;
Te2 = 60;  
Ti = 2000;
S = 0.4;
% Ions
omega_pi = irf_plasma_calc(B,n,no,Te1,Ti,'Fpp')*2*pi; % rad/s
vthi = irf_plasma_calc(B,n,no,Te1,Ti,'Vtp'); % m/s 
vdi = 0;
% Background electrons
omega_pe1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Fpe')*2*pi; % rad/s
vthe1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Vte'); % m/s 
vde1 = 0; % m/s
% Beam electrons
omega_pe2 = irf_plasma_calc(B,n2,no,Te2,Ti,'Fpe')*2*pi; % rad/s
vthe2 = irf_plasma_calc(B,n2,no,Te2,Ti,'Vte'); % m/s 
vde2 = S*vthe1;
% Other
lamD = irf_plasma_calc(B,n,no,Te1,Ti,'Ld'); % m
omega_bune = omega_pe1^(1/3)*omega_pi^(2/3)*16^(-1/3);

% Make dimensionless parameters
om_norm = omega_pe1;
om_i = omega_pi/om_norm;
om_e1 = omega_pe1/om_norm;
om_e2 = omega_pe2/om_norm;

vt_norm = vthe1;
vt_i = vthi/vt_norm;
vt_e1 = vthe1/vt_norm;
vt_e2 = vthe2/vt_norm;
vd_e2 = vde2/vt_norm;

disp(['o_bg    = ',num2str(om_e1)])
disp(['o_beam  = ',num2str(om_e2)])
disp(['o_i     = ',num2str(om_i)])
disp(['vt_bg   = ',num2str(vt_e1)])
disp(['vt_beam = ',num2str(vt_e2)])
disp(['vt_i    = ',num2str(vt_i)])
disp(['vd_e2   = ',num2str(vd_e2)])
