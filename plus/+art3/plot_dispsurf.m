%% Figure for lic-presentation
% Good from above, #5: z - omega, c - gamma, and #7: z - vph, c - gamma.
figure(28)

[gS,gk]=meshgrid(S,k);
plot_imag = x_imag_store;
Slim = 0.81;
klim = 0.11;
plot_imag(gS>Slim) = NaN;
plot_imag(gk<klim) = NaN;
plot_imag(x_real_store<0) = NaN;
plot_real = x_real_store;
plot_real(gS>Slim) = NaN;
plot_real(gk<klim) = NaN;
plot_real(x_real_store<0) = NaN;


kmat = repmat(k,nv,1)';
plot_v = plot_real./kmat;

omega_pe = omega_pe1;

cmap = cn.cmap('islands');
colormap(cmap)
nPanels = 2;
for ii=1:nPanels; h(ii)=subplot(1,ceil(nPanels/1),ii); end
set(gcf,'position',[268   579   842   257]);
isub=0;
axpos = cell(nPanels,1);

if 1
    isub=isub+1; hca = h(isub); 
    toplot = plot_imag'/omega_pe;
    toplot(toplot<-0.04) = NaN;
    surf(hca,k,S,toplot)
    axpos{isub} = get(hca,'position');
    xlabel(hca,'k\lambda_D')
    ylabel(hca,'v_d/v_{te1}')
    zlabel(hca,'\omega_i/\omega_{pe}')
    %hc(isub) = colorbar('peer',hca);
    set(hca,'clim',(max(get(hca,'zlim')))*[-1 1],'zlim',0.04*[-1 1])
end
if 0
    isub=isub+1; hca = h(isub); 
    surf(hca,k,S,plot_real'/omega_pe,plot_imag'/omega_pe)
    axpos{isub} = get(hca,'position');
    xlabel(hca,'k\lambda_D')
    ylabel(hca,'v_d/v_{te1}')
    zlabel(hca,'\omega_{r}/\omega_{pe}')
    %hc(isub) = colorbar('peer',hca);
    %ylabel(hc(isub),'f_i [Hz]')
    %set(hc(isub),'ylim',max(get(hc(isub),'ylim'))*[-1 1])
    set(hca,'clim',max(get(hca,'clim'))*[-1 1])
    %set(hca,'zaxis',[0 2000])
end
if 1
    isub=isub+1; hca = h(isub); 
    surf(hca,k,S,(plot_v)'*lamD/1000,plot_imag'/omega_pe)    
    axpos{isub} = get(hca,'position');
    xlabel(hca,'k\lambda_D')
    ylabel(hca,'v_d/v_{te1}')
    zlabel(hca,'v_{ph} [km/s]')
    hc(isub) = colorbar('peer',hca);
   
    ylabel(hc(isub),'\omega_i/\omega_{pe}');
    %set(hca,'zscale','log');
    set(hca,'clim',max(get(hca,'clim'))*[-1 1],'zlim',[0 12000])
    %set(hca,'clim',max(get(hc(isub),'ylim'))*[-1 1])
    %set(hca,'zaxis',[0 2000])
end


displace = {[0 0 -0.05 0],[-0.05 0 -0.05 0],[0.05 0 0 0]};
for kk = 1:2%nPanels
    %shading(h(kk),'flat')
    set(h(kk),'xlim',[max([klim k(1)]) k(end)],'ylim',[0.2 Slim],...
        'clim',max(get(h(2),'clim'))*[-1 1],...
        'position',axpos{kk}+displace{kk})
end
%set(h(1),'zlim',[-0.08 0.07]);
%% 
figure(29)

f = @(n,v,vd,vt) (n^(1)*1e6/((pi^(1/2)*vt)^(3))).*exp(-(v-vd).^2/(1*vt.^2));
units = irf_units;    
v = linspace(-2,2,500)*vthe1;
fiscale = 0.0004;
plot_vde2 = [1 0.8];
%pind(1) = find(S>plot_vde2(1),1,'first');
%pind(2) = find(S>plot_vde2(2),1,'first');
plot(v/vthe1,f(n1,v,vde1,vthe1)+f(n2,v,vthe1*plot_vde2(1),vthe2),...
    v/vthe1,f(n1,v,vde1,vthe1)+f(n2,v,vthe1*plot_vde2(2),vthe2),...
    v/vthe1,fiscale*f(n,v,vdi,vthi),...
    'linewidth',2);
axpos{isub} = get(hca,'position');
%t('f(v_{||})')
xlabel(hca,'v_{||}/v_{te}')
%irf_legend({' ','f_i'},[0.38 0.8])
%irf_legend({'f_e',' '},[0.7 0.5])
%set(gca,'yscale','log','xscale','log','ylim',[1e-30 1e-15])

