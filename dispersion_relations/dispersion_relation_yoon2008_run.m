% dispersion_relation_yoon2008_run.m

% prepare input
%tint = toepoch([2008 04 22 17 37 01.2;2008 04 22 17 37 02.0])';
%tint = toepoch([2008 04 22 17 37 13.0;2008 04 22 17 37 14.0])'; % em waves too in high beta
%tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])'; % em waves too in high beta
%tool.run_single_event
plot_option = '2D';

units = irf_units;
kB = 1.381e-23; % J/K
e = 1.602e-19; % C

%Te_loc_MK = irf_resamp(parTe1MK,mean(tint),'nearest'); Te_loc_MK = Te_loc_MK(2);
%Ti_loc_MK = irf_resamp(Ti1MK,mean(tint),'nearest'); Ti_loc_MK = Ti_loc_MK(2);
%Te_loc = Te_loc_MK*units.kB/units.e*1e6;
%Ti_loc = Ti_loc_MK*units.kB/units.e*1e6;


guide_field = 0.00001; % fraction of B0
B0 = B0;
n = n_loc;
no = 0;
Te = Te_loc; TeK = units.e/units.kB*Te;
Ti = Ti_loc;  TiK = units.e/units.kB*Ti;
LnRop = 0.5;

% Calculate all the needed parameters
omega_pe = irf_plasma_calc(B0,n,no,Te,Ti,'Fpe')*2*pi; % rad/s
omega_ce = irf_plasma_calc(B0,n,no,Te,Ti,'Fce')*2*pi; % rad/s
omega_pi = irf_plasma_calc(B0,n,no,Te,Ti,'Fpp')*2*pi; % rad/s
omega_ci = irf_plasma_calc(B0,n,no,Te,Ti,'Fcp')*2*pi; % rad/s
omega_lh = irf_plasma_calc(B0,n,no,Te,Ti,'Flh')*2*pi; % rad/s
vthi = irf_plasma_calc(B0,n,no,Te,Ti,'Vtp')/sqrt(2); % m/s sqrt(kBTe/me)
vthe = irf_plasma_calc(B0,n,no,Te,Ti,'Vte')/sqrt(2); % m/s sqrt(kBTe/me) 
vthi = irf_plasma_calc(B0,n,no,Te,Ti,'Vtp'); % m/s sqrt(2kBTe/me)
vthe = irf_plasma_calc(B0,n,no,Te,Ti,'Vte'); % m/s sqrt(2kBTe/me) 
Rop = irf_plasma_calc(B0,n,no,Te,Ti,'Rop'); % m
Ln = LnRop*Rop;                    % gradient length scale, m
vdi = kB*TiK/e/(B*1e-9)/Ln;     % ion diamagnetic drift, m/s
vde = -kB*TeK/e/(B*1e-9)/Ln;    % electron diamagnetic drift, m/s

disp(['vdi/vti=' num2str(vdi/vthi,'%.1f') ', Ti/Te=' num2str(Ti/Te,'%.1f')])

% Put all the physical parameters into the dispersion function
fv = [omega_pi,omega_pe,omega_ci,omega_ce,omega_lh,vthe,vthi,vdi,vde,Rop,Ln,guide_field];

% Set up figure for drawing results
switch plot_option
    case 'noplot'
    otherwise
        figure(22)
        for nPlot = 1:2     
            h(nPlot) = subplot(1,2,nPlot); hold(h(nPlot),'on');        
            xlabel(h(nPlot),'k_z  v_{th,e}/ \omega_{ce}','Fontsize',14);    
            set(h(nPlot),'XGrid','on','YGrid','on','Fontsize',14);
        end
        ylabel(h(1),'\omega_r/\omega_{LH}','Fontsize',14);
        ylabel(h(2),'\omega_i/\omega_{LH}','Fontsize',14);
end
% Initial guess
x_real_store2=0.03*omega_lh % good for myfun_20070831
x_imag_store2=0.01*omega_lh

% k*rho_e to loop over
kx = 0;0.1:0.1:1;   % loop over kx*rho_e.
ky = 0.1:0.05:2;   % loop over ky*rho_e.
kz = 0;   % loop over kz*rho_e

for ix = 1:numel(kx)
    for iy = 1:numel(ky)
        for iz = 1:numel(kz)
            % Do solver
            af = @(temp) dispersion_relation_yoon2008(temp,kx(ix),ky(iy),kz(iz),fv);
            [x,FVAL,EXITFLAG] = fsolve(af,[x_real_store2(end,end) x_imag_store2(end,end)], optimset('GradObj','on','display','off'));        

            % Store new results in array
            x_real_store2(ix,iy,iz)=x(1);
            x_imag_store2(ix,iy,iz)=x(2);

            % Draw results
            switch plot_option
                case 'noplot'
                case '2D'
                    plot(h(1),real(ky(iy))/sqrt(2),x(1)/omega_lh,'r+');
                    plot(h(2),real(ky(iy))/sqrt(2),x(2)/omega_lh,'r+');
                    set(h(1),'ylim',[0 4])
                    set(h(2),'ylim',[-1 1])
                    drawnow
                case '3D'
                    plot3(h(1),real(kx(ix)),real(ky(iy)),x(1)/omega_lh,'r+'); 
                    plot3(h(2),real(kx(ix)),real(ky(iy)),x(2)/omega_lh,'r+');
                    set(h(1),'zlim',[0 5])
                    set(h(2),'zlim',[0 5])
                    drawnow
            end            
        end
    end
end                                                       