%% Prepare xes1 input file
Tebg = 2000; 
Teb = 60;
Ti = 2000; 

B = 25;
n = 0.06;
no = 0;

% Ions
omega_pi = irf_plasma_calc(B,n,no,Tebg,Ti,'Fpp')*2*pi; % rad/s
vthi = irf_plasma_calc(B,n,no,Tebg,Ti,'Vtp'); % m/s 
vdi = 0;
% Total electron plasma frequency
omega_pe = irf_plasma_calc(B,n,no,Tebg,Ti,'Fpe')*2*pi; % rad/s


units = irf_units;
%RR = [0.1 0.15 0.2 0.2 0.2 0.2 0.4 0.5  0.5  0.5 0.6 0.98]; 
%SS = [0.4 0.7 0.7 0.4 1.5 0.5 0.5 0.55 0.65 0.8 0.6 0.9 ];
RR = [0.3 0.4 0.5]; 
SS = [0.45 0.5 0.8];
disp('---- iPIC3D input ----')
v_norm = units.c;

for ir = 1:numel(RR)    
    vte1 = irf_plasma_calc(B,n,no,Tebg,Ti,'Vte'); % m/s 
    vte2 = irf_plasma_calc(B,n,no,Teb,Ti,'Vte'); % m/s 
    vti = irf_plasma_calc(B,n,no,Tebg,Ti,'Vtp'); % m/s 
    vde2 = SS(ir)*vte1;

    iPIC_vte1 = vte1/v_norm;
    iPIC_vte2 = vte2/v_norm;
    iPIC_vti = vti/v_norm;
    iPIC_vde1 = 0/v_norm;
    iPIC_vde2 = vde2/v_norm;  
    iPIC_vdi = 0/v_norm;
    
    disp(['ne1 = ' num2str((1-RR(ir)),'%.2f') ', ne2 = ' num2str((RR(ir)),'%.2f') ', ni = 1.00, ',...
          'vte1 = ' num2str(iPIC_vte1,'%.4f') ', vte2 = ' num2str(iPIC_vte2,'%.4f') ', vti = ' num2str(iPIC_vti,'%.4f') ', ',...
          'vde1 = ' num2str(iPIC_vde1,'%.4f') ', vde2 = ' num2str(iPIC_vde2,'%.4f') ', vdi = ' num2str(iPIC_vdi,'%.4f')])   
end