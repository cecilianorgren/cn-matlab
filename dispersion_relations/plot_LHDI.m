% Run first solver_LHDI
ts = 16;
set(gcf,'defaultAxesFontSize',ts);
set(gcf,'defaultTextFontSize',ts);
set(gcf,'position',[560 647 490 301]);

flh = omega_lh/2/pi;

plot(aaa,x_real_store2/omega_lh,aaa,x_imag_store2/omega_lh,'linewidth',2)
ylabel('\omega/\omega_{LH}','Fontsize',ts)
%xlabel('k_z  v_{th,e}/ \omega_{ce}','Fontsize',ts);
xlabel('k_y  \rho_e','Fontsize',ts);
%legend('\omega_r','\omega_i')

ind=4;
text(aaa(end-ind),1.1*x_real_store2(end-ind)/omega_lh,'\omega_{r}','Fontsize',ts)
text(aaa(end-ind),1.2*x_imag_store2(end-ind)/omega_lh,'\omega_{i}','Fontsize',ts)
text(0.1,0.9*max(x_real_store2/omega_lh),['f_{LH} = ' num2str(flh,'%.0f') ' Hz'],'Fontsize',ts)

%set(gca,'ylim',[0 max(get(gca,'ylim'))])
set(gca,'ylim',[0 1])