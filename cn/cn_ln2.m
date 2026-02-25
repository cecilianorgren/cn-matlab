function ln = cn_ln2(t1,t2,t1str,t2str,gsmB3,gsmB4,Te,Ti,Ne,Ni,dx,dy,dz,M,r_i)
mu0=4*pi*1e-7;
kb=1.38e-23;                % J/K
T=(Te+Ti)*1e6/(8.61734*10)     % K
n=(Ne+Ni)*1e6 /2;                   % m^-3

c_eval('gseB?=c_coord_trans(''gsm'',''gse'',gsmB?,''CL_ID'',?);',3:4);

dB=irf_add(1,gseB3,-1,gseB4);
dB=cn_toepoch(t1,t2,dB);
gseB3=cn_toepoch(t1,t2,gseB3);
gseB4=cn_toepoch(t1,t2,gseB4);

prel=M*dB(:,2:4)';
bdB=[dB(:,1) prel'];
prel=M*gseB3(:,2:4)';
bB3=[gseB3(:,1) prel'];
prel=M*gseB4(:,2:4)';
bB4=[gseB4(:,1) prel'];

%bB4b=cn_xyz(gseB4(:,1:4),[gseB4(1,1) M(3,:)],[gseB4(1,1) M(2,:)]);

bB3=irf_resamp(bB3,bdB);
bB4=irf_resamp(bB4,bdB);
%figure;irf_plot({bB3,bB4},'comp')
bdB=irf_abs(bdB);
j=cn_j(bdB,[-dx -dy -dz]);
jk=[bdB(:,1) -bdB(:,4)*1e-9/(mu0*dy*1e3)]; % A/m^2
jnorm=[bdB(:,1) bdB(:,4)*1e-9/(mu0*dx*1e3)];
jxB=irf_multiply(1,jk,1,bB3(:,[1 4]),1);
meanbB=[bB3(:,1) mean([bB3(:,4) bB4(:,4)],2)];
if 1
    figure;irf_plot(bdB)
end
if 1 % plot current
    figure('name','Current from C3-C4');
    irf_plot(jk);hold on;
    irf_plot(jnorm,'g');hold on;
    legend('j_k','j_{norm}')
    title('Current from C3-C4')
    %ln=[jxB(:,1) (0.1e-6*30)./jxB(:,2)];
    %ln=[jxB(:,1) (0.1e-6*10e6*kb)./jxB(:,2)];
end
lnoverr_i=[j(:,1) abs(kb*T*n./(j(:,2).*meanbB(:,2)*1e-9*r_i*1000))]; 
if 0 % plot ln over ri
figure('name','Ln');%irf_plot(jk);hold on;
%irf_plot(jnorm,'g');hold on;
irf_plot(lnoverr_i,'g');hold on;
title(['Gradient length scale L_n  ( r_i = ' num2str(r_i,'%.0f') ' km )'])
ylabel('L_n/\rho_i')
set(gcf,'PaperPositionMode','auto');
eval(['print -depsc2 ',t1str,'_',t2str,'_ln.eps']);
end

ln=lnoverr_i;