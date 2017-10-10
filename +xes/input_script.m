Te1 = 2000;
Te2 = 60;
Ti = 2000;
n = 1;
no = 0;
B=10;

omega_pe = irf_plasma_calc(B,n,no,Te1,Ti,'Fpe')*2*pi; % m/s 
%% Prepare xes1 input file
units = irf_units;
RR = [0.1 0.15 0.2 0.2 0.2 0.2 0.4 0.5  0.5  0.5 0.6 0.98]; 
SS = [0.4 0.7 0.7 0.4 1.5 0.5 0.5 0.55 0.65 0.8 0.6 0.9 ];
disp('----')
for ir = 1:numel(RR)    
    ope1 = omega_pe*sqrt(1-RR(ir));
    ope2 = omega_pe*sqrt(RR(ir));
    opi = omega_pe*sqrt(units.me/units.mp);
    vte1 = irf_plasma_calc(B,n,no,Te1,Ti,'Vte'); % m/s 
    vte2 = irf_plasma_calc(B,n,no,Te2,Ti,'Vte'); % m/s 
    vti = irf_plasma_calc(B,n,no,Te1,Ti,'Vtp'); % m/s 
    vde2 = SS(ir)*vte1;

    xes1_ope1 = ope1/omega_pe;
    xes1_ope2 = ope2/omega_pe;
    xes1_opi = opi/omega_pe;
    xes1_vte1 = vte1/vte1;
    xes1_vte2 = vte2/vte1;
    xes1_vti = vti/vte1;
    xes1_vde2 = vde2/vte1;
    
    disp(['R = ' num2str(RR(ir),'%.2f') ', S = ' num2str(SS(ir),'%.2f') ', wpe1 = ' num2str(xes1_ope1,'%.3f') ', wpe2 = ' num2str(xes1_ope2,'%.3f') ', wpi = ' num2str(xes1_opi,'%.3f') ', vte1 = ' num2str(xes1_vte1,'%.3f') ', vte2 = ' num2str(xes1_vte2,'%.3f') ', vti = ' num2str(xes1_vti,'%.3f') ', vde2 = ' num2str(xes1_vde2,'%.3f')])
end
