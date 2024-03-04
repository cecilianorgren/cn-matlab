%% Draw distribution with trapped passing boundary
units = irf_units;

% Boundary conditions
% Parallel potential
phi = 20; % eV
% Magnetic field
B1 = 10e-9; % T 
B2 = 3e-9; % T 

% Upstream distibution
T1 = 50; % eV

%Ek_perp = 50; 
%vperp1 = sqrt(Ek_perp*units.eV*2/units.me);


% Ek = Ekpar + Ekperp;

% Without phi
% mu*B1 = mu*B2

vmax = sqrt(200*units.eV*2/units.me);


vpar = linspace(-vmax,vmax,200);
vperp = linspace(0,vmax,100);

[VPAR,VPERP] = ndgrid(vpar,vperp);
VABS = sqrt(VPAR.^2 + VPERP.^2);
EK2 = (VPAR.^2 + VPERP.^2)*units.me/2;
EKpar2 = (VPAR.^2)*units.me/2;

MU2 = units.me*VPERP.^2/(2*B2);
MU1 = MU2;

EKpar1 = EK2 - MU1*B1 - phi*units.e;
EKperp1 = MU1*B1;

EKpar1_lim1 = EKpar2 - 0*MU1*B1 - phi*units.e;
EKpar1_lim2 = EK2 - MU1*B1 - 0*phi*units.e;
EKperp1 = MU1*B1;

VPAR1 = sqrt(2*EKpar1/units.me);
VPERP1 = sqrt(2*EKperp1/units.me);


vt1 = sqrt(T1*units.eV*2/units.me);
f3 = @(v,vt,n,vd) n/((pi)^(3/2)*vt.^3)*exp(-(v-vd).^2./vt./vt);
f2 = @(v,vt,n,vd) n/((pi)^(2/2)*vt.^2)*exp(-(v-vd).^2./vt./vt);
f1 = @(v,vt,n,vd) n/((pi)^(1/2)*vt.^1)*exp(-(v-vd).^2./vt./vt);

f2_ = @(v1,v2,vt,n,vd) n/((pi)^(2/2)*vt.^2)*exp(-(v1-vd).^2./vt./vt-(v2-vd).^2./vt./vt);

F1 = f2(VABS,vt1,1,0);
F2 = VABS*0;
ipass = EKpar1>0;
itrap = EKpar1<0;
%ipass = EKpar1>-1e20;
%itrap = EKpar1<-1e20;
F2(ipass) = f2_(VPAR1(ipass),VPERP1(ipass),vt1,1,0);
F2(itrap) = f2_(0,VPERP1(itrap),vt1,1,0);

h = setup_subplots(3,2,'horizontal');
isub = 1;

hca = h(isub); isub = isub + 1;
contourf(hca,VPAR*1e-3,VPERP*1e-3,EK2/units.eV)
shading(hca,'flat')
hb = colorbar(hca);
hb.YLabel.String = 'E_{k}^2';


hca = h(isub); isub = isub + 1;
contourf(hca,VPAR*1e-3,VPERP*1e-3,EKpar1/units.eV,15)
shading(hca,'flat')
hb = colorbar(hca);
colormap(hca,irf_colormap('waterfall'))
hca.CLim = max(abs(EKpar1(:))/units.eV)*[-1 1];
hold(hca,'on')
contour(hca,VPAR*1e-3,VPERP*1e-3,EKpar1/units.eV,[0 0],'k','LineWidth',1)
hold(hca,'off')
hb.YLabel.String = 'E_{k||}^1';

hca = h(isub); isub = isub + 1;
contourf(hca,VPAR*1e-3,VPERP*1e-3,F1,15)
shading(hca,'flat')
hb = colorbar(hca);
hb.YLabel.String = 'f^1';
hca.XLabel.String = 'v_{||} (km/s)';
hca.YLabel.String = 'v_{\perp} (km/s)';
%linkprop(h(isub+[-1 -2]),{'CLim'})
irf_legend(hca,{['B_1 = ' num2str(B1*1e9) ' nT'],['\phi_1 = 0 eV']}',[0.02 0.98],'color','k','fontsize',12)
axis(hca,'equal')

hca = h(isub); isub = isub + 1;
%pcolor(hca,VPAR*1e-3,VPERP*1e-3,F2)
contourf(hca,VPAR*1e-3,VPERP*1e-3,F2,15)
shading(hca,'flat')
hb = colorbar(hca);
hold(hca,'on')
contour(hca,VPAR*1e-3,VPERP*1e-3,EKpar1/units.eV,[0 0],'k','LineWidth',1)
contour(hca,VPAR*1e-3,VPERP*1e-3,EKpar1_lim1/units.eV,[0 0],'k:','LineWidth',1)
contour(hca,VPAR*1e-3,VPERP*1e-3,EKpar1_lim2/units.eV,[0 0],'k:','LineWidth',1)
hold(hca,'off')
hb.YLabel.String = 'f^2';
hca.XLabel.String = 'v_{||} (km/s)';
hca.YLabel.String = 'v_{\perp} (km/s)';
irf_legend(hca,{['B_2 = ' num2str(B2*1e9) ' nT'],['\phi_2 = ' num2str(phi) ' eV']}',[0.02 0.98],'color','k','fontsize',12)
axis(hca,'equal')
%irf_legend(hca,{['\mathcal{E} = ' num2str(B2*1e9) ' nT'],['\phi_2 = ' num2str(phi) ' eV']}',[0.02 0.98],'color','k','fontsize',12)
irf_legend(hca,{'E_{||,1}= E_{||,2}- e\phi' }',[0.99 0.98],'color','k','fontsize',12)
irf_legend(hca,{'E_{||,1}= E_{||,2}- \mu B_1' }',[0.99 0.8],'color','k','fontsize',12)

hca = h(isub); isub = isub + 1;
pcolor(hca,VPAR*1e-3,VPERP*1e-3,real(VPAR1))
shading(hca,'flat')
hb = colorbar(hca);
hb.YLabel.String = 'v_{||}^1';

hca = h(isub); isub = isub + 1;
pcolor(hca,VPAR*1e-3,VPERP*1e-3,real(VPERP1))
shading(hca,'flat')
hb = colorbar(hca);
hb.YLabel.String = 'v_{\perp}^1';


%colormap(irf_colormap('candy4'))
colormap(pic_colors('waterfall'))
colormap(pic_colors('candy4'))