% Script for solving disersion relation for the LHDI.
% The physical variables are assigned in myfun_20070831.
% x: normal direction, y: drift direction, z: magentic field direction
global omega_lh

% set up figure for drawing results
figure(2)    
for nPlot = 1:2     
    h(nPlot) = subplot(1,2,nPlot); hold(h(nPlot),'on');        
    xlabel(h(nPlot),'k_z  v_{th,e}/ \omega_{ce}','Fontsize',14);    
    set(h(nPlot),'XGrid','on','YGrid','on','Fontsize',14);
end
ylabel(h(1),'\omega_r (Hz)','Fontsize',14);
ylabel(h(2),'\omega_i (Hz)','Fontsize',14);

% initial guess
omega_lh=300;
x_real_store2=0.043*omega_lh*0.1; % good for myfun_20070831
x_imag_store2=0.015*omega_lh*0.1;

% k*rho_e to loop over
ky = 0.01:0.05:2;   % outer loop, over ky*rho_e.
kz = [0];          % inner loop over kz*rho_e

for iy = 1:numel(ky)
    for iz = 1:numel(kz)
        af = @(temp) myfun_20070831(temp,0,ky(iy),kz(iz));               
        [x,FVAL,EXITFLAG] = fsolve(af,[x_real_store2(end,end) x_imag_store2(end,end)], optimset('GradObj','on','display','off'));        
    
        % store new results in array 
        x_real_store2(iy,iz)=x(1);
        x_imag_store2(iy,iz)=x(2);
                       
        % draw results        
        plot(h(1),real(ky(iy)),x(1)/(2*pi),'r+'); 
        plot(h(2),real(ky(iy)),x(2)/(2*pi),'r+');                                                
        drawnow                       
    end
end