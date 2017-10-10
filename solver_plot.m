%% surf plot
figure(17)
plot_imag = x_imag_store;
plot_imag(x_real_store<0) = NaN;
plot_real = x_real_store;
plot_real(x_real_store<0) = NaN;

kmat = repmat(k,nv,1)';

plot_v = plot_real./kmat;

omega_pe = omega_pe1;

nPanels = 4;
for ii=1:nPanels; h(ii)=subplot(2,2,ii); end

hca = h(1);
surf(hca,k,S,plot_real'/omega_pe)
xlabel(hca,'k\lambda_D')
ylabel(hca,'v_d/v_{te1}')
zlabel(hca,'\omega_r/\omega_{pe}')

hca = h(2);
surf(hca,k,S,plot_imag'/omega_pe)
xlabel(hca,'k\lambda_D')
ylabel(hca,'v_d/v_{te1}')
zlabel(hca,'\omega_i/\omega_{pe}')

hca = h(3);
surf(hca,k,S,plot_v'*lamD/1000)
xlabel(hca,'k\lambda_D')
ylabel(hca,'v_d/v_{te1}')
zlabel(hca,'v_{ph}')
%set(hca,'zaxis',[0 2000])

hca = h(4);
surf(hca,k,S,plot_imag'/omega_pe,plot_v'*lamD/1000)
xlabel(hca,'k\lambda_D')
ylabel(hca,'v_d/v_{te1}')
zlabel(hca,'v_{ph} [km/s]')
hc(4) = colorbar('peer',hca);
%set(hca,'clim',[0 3500])

%% pcolor plot

omega_pe=omega_pe1;

if 0 % clean up results
    toplim = omega_pe*0.9*2*pi;
    plot_imag = x_imag_store;
    plot_imag(x_real_store<0) = NaN; % take away negative frequencies
    plot_imag(x_real_store>toplim) = NaN;
    
    plot_real = x_real_store;
    plot_real(x_real_store<0) = NaN;
    plot_real(x_real_store>toplim) = NaN;
    kmat = repmat(k,nv,1)';
    plot_v = plot_real./kmat;
end


if 1
    normalization = 2*pi;
    f_string = ' [s^{-1}]';
    flim = [0 0.5*omega_pe/2/pi];
else
    normalization = omega_pe;
    f_string = '\omega_r/\omega_{pe}';
    flim = [0 1];
end

figure(20)
set(gcf,'position',[113   104   444   852])
nPanels = 3;
for ii=1:nPanels; h(ii)=subplot(3,1,ii); end

hca = h(1);
pcolor(hca,k,S,plot_real'/normalization)
hc(1) = colorbar('peer',hca);
ylabel(hc(1),['f_r  ' f_string])
%set(hca,'clim',flim)

hca = h(2);
pcolor(hca,k,S,plot_imag'/normalization)
hc(2) = colorbar('peer',hca);
ylabel(hc(2),['f_i  ' f_string])
%set(hca,'clim',max(get(hca,'clim'))*[-1 1]*1)
%set(hca,'clim',[-1 1]*20)
%set(hc(2),'ytick',-100:2:100)
%cmp = irf_colormap('poynting');
%colormap(cmp)

hca = h(3);
pcolor(hca,k,S,(plot_v)'*lamD/1000)
hc(3) = colorbar('peer',hca);
ylabel(hc(3),'v_{ph}  [km/s]')
%set(h(3),'clim',[0 1500])

for ii=1:3; 
    shading(h(ii),'flat'); 
    xlabel(h(ii),'k\lambda_D')
    ylabel(h(ii),'v_d/v_{te1}')
    hold(h(ii),'off')
    axis(h(ii),'square')
end
colormap(cn.cmap('islands'))
title(h(1),['f_{pe} = ' num2str(omega_pe1/2/pi,'%.0f') ' Hz,  v_{te1} = ' num2str(vthe1*1e-6,'%.1f') '\cdot10^3 km/s,  R = ' num2str(R,'%.2f')])

%% Different presentations of dispersion relation
figure(28)

[gS,gk]=meshgrid(S,k);
plot_imag = x_imag_store;
Slim = 1.8;0.85;
klim = 0.12;
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
nPanels = 7;
for ii=1:nPanels; h(ii)=subplot(2,ceil(nPanels/2),ii); end
isub=0;
if 1
    isub=isub+1; hca = h(isub); 
    surf(hca,k,S,plot_real'/omega_pe)
    xlabel(hca,'k\lambda_D')
    ylabel(hca,'v_d/v_{te1}')
    zlabel(hca,'\omega_r/\omega_{pe}')
end
if 1
    isub=isub+1; hca = h(isub); 
    surf(hca,k,S,plot_imag'/omega_pe)
    xlabel(hca,'k\lambda_D')
    ylabel(hca,'v_d/v_{te1}')
    zlabel(hca,'\omega_i/\omega_{pe}')
    hc(isub) = colorbar('peer',hca);
    set(hca,'clim',(max(get(hca,'zlim')))*[-1 1])
end
if 1
    isub=isub+1; hca = h(isub); 
    surf(hca,k,S,(plot_v)'*lamD/1000)
    xlabel(hca,'k\lambda_D')
    ylabel(hca,'v_d/v_{te1}')
    zlabel(hca,'v_{ph} [km/s]')
    hc(isub) = colorbar('peer',hca);
    %set(hc(isub),'yscale','log');
    set(hca,'zscale','log');
    %set(hca,'zaxis',[0 2000])
end
if 1
    isub=isub+1; hca = h(isub); 
    surf(hca,k,S,plot_imag'/omega_pe,plot_v'*lamD/1000)
    xlabel(hca,'k\lambda_D')
    ylabel(hca,'v_d/v_{te1}')
    zlabel(hca,'\omega_i/\omega_{pe}')
    hc(isub) = colorbar('peer',hca);
    ylabel(hc(isub),'v_{ph} [km/s]')
    %set(hca,'clim',[0 3500])
end
if 1
    isub=isub+1; hca = h(isub); 
    surf(hca,k,S,plot_real'/omega_pe,plot_imag'/2/pi)
    xlabel(hca,'k\lambda_D')
    ylabel(hca,'v_d/v_{te1}')
    zlabel(hca,'\omega_{r}/\omega_{pe}')
    hc(isub) = colorbar('peer',hca);
    ylabel(hc(isub),'f_i [Hz]')
    set(hc(isub),'ylim',max(get(hc(isub),'ylim'))*[-1 1])
    set(hca,'clim',max(get(hc(isub),'ylim'))*[-1 1])
    %set(hca,'zaxis',[0 2000])
end
if 1
    isub=isub+1; hca = h(isub); 
    surf(hca,k,S,plot_imag'/omega_pe,plot_real'/2/pi)
    xlabel(hca,'k\lambda_D')
    ylabel(hca,'v_d/v_{te1}')
    zlabel(hca,'\omega_i/\omega_{pe}')
    hc(isub) = colorbar('peer',hca);
    ylabel(hc(isub),'f_r [Hz]')
    %set(hca,'clim',[0 3500])
end
if 1
    isub=isub+1; hca = h(isub); 
    surf(hca,k,S,(plot_v)'*lamD/1000,plot_imag'/omega_pe)
    xlabel(hca,'k\lambda_D')
    ylabel(hca,'v_d/v_{te1}')
    zlabel(hca,'v_{ph} [km/s]')
    hc(isub) = colorbar('peer',hca);
   
    ylabel(hc(isub),'\omega_i/\omega_{pe}');
    %set(hca,'zscale','log');
    set(hca,'clim',max(get(hc(isub),'ylim'))*[-1 1])
    %set(hca,'zaxis',[0 2000])
end

for kk = 1:nPanels
    %shading(h(kk),'flat')
    set(h(kk),'xlim',[klim k(end)],'ylim',[S(1) min([Slim S(end)])])
end

%% Figure for lic-presentation
% Good from above, #5: z - omega, c - gamma, and #7: z - vph, c - gamma.
figure(28)

[gS,gk]=meshgrid(S,k);
plot_imag = x_imag_store;
Slim = 0.85;
klim = 0.22;
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
nPanels = 3;
for ii=1:nPanels; h(ii)=subplot(1,ceil(nPanels/1),ii); end
isub=0;

axpos = cell(nPanels,1);

if 1
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
    set(hca,'clim',max(get(hca,'clim'))*[-1 1])
    %set(hca,'clim',max(get(hc(isub),'ylim'))*[-1 1])
    %set(hca,'zaxis',[0 2000])
end
if 0 % Add distribution functions
    isub=isub+1; hca = h(isub); 
    f = @(n,v,vd,vt) (n^(1)/((pi^(1/2)*vt)^(3))).*exp(-(v-vd).^2/(1*vt.^2));
    units = irf_units;    
    v = linspace(-5,5,500)*vthe1;
    fiscale = 0.01;
    plot_vde2 = [0.5 0.7];
    pind(1) = find(S>plot_vde2(1),1,'first');
    pind(2) = find(S>plot_vde2(2),1,'first');
    plot(hca,v/vthe1,f(n1,v,vde1,vthe1)+f(n2,v,vde2(pind(1)),vthe2),...
        v/vthe1,f(n1,v,vde1,vthe1)+f(n2,v,vde2(pind(2)),vthe2),...
        v/vthe1,fiscale*f(n,v,vdi,vthi),...
        'linewidth',2);
    %t('f(v_{||})')
    xlabel(hca,'v_{||}/v_{te}')
    %irf_legend({' ','f_i'},[0.38 0.8])
    %irf_legend({'f_e',' '},[0.7 0.5])
end

for kk = 1:nPanels
    %shading(h(kk),'flat')
    set(h(kk),'xlim',[klim k(end)],'ylim',[S(1) Slim],...
        'position',axpos{kk})
end
