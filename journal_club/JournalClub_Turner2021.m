B_nT = linspace(10,50000,50);
B_nT = logspace(log10(10),log10(5000),50);
B_G = B_nT*1e-9*1e4;

Emax = 0.6; % MeV
hca = subplot(1,1,1);
loglog(hca,B_nT,Emax*B_G.^(-1))
ylim = hca.YLim;
hold(hca,'on')
Muall = [1764,879,439,219,109];
for iMu = 1:numel(Muall)
  B_tmp = Emax/Muall(iMu)*1e5;
  %loglog(hca,B_nT,Muall(iMu)+B_nT*0,'k--')
  loglog(hca,B_tmp+[0 0],hca.YLim,'k--')
  loglog(hca,B_tmp,Muall(iMu),'*')
  %loglog(hca,B_nT,879+B_nT*0,'k--')
  %loglog(hca,B_nT,439+B_nT*0,'k--')
  %loglog(hca,B_nT,219+B_nT*0,'k--')
  %loglog(hca,B_nT,109+B_nT*0,'k--')
  %text(hca,B_nT(end)*0.8,Muall(iMu),sprintf('Mu = %4.0f MeV/G',Muall(iMu)),'horizontalalignment','right','verticalalignment','bottom','fontsize',12)
  text(hca,B_tmp,Muall(iMu),sprintf('  Mu = %4.0f MeV/G',Muall(iMu)),'horizontalalignment','left','verticalalignment','bottom','fontsize',12)
end
hold(hca,'off')

hca.YLim = ylim;
hca.XLim = B_nT([1 end]);
hca.FontSize = 14;


xlabel('B (nT)')
ylabel('Mu (MeV/G)')
title('Mu cut-off: Mu_{max} = E_{max}/B')

%c_eval('feeps?e = mms.get_data(''Omnifluxelectron_epd_feeps_brst_l2'',tint,?);',ic) 
%% Energy vs velocity
% E = m*c^2 - both rest mass and kinetic energy
% Kinetic energy EK = m*c^2 - m0*c^2
% v = sqrt(2*K/m0)
units = irf_units;
E = logspace(-1,log10(2.99e6),90); % eV
%ER_SI = units.me*units.c^2;
%ER_eV = ER_SI/units.eV;
%EK = E - ER_eV;
%v = sqrt(EK*units.eV*2/units.me); % m/s

v = units.c*sqrt(1-1./(E.*units.e./(units.me*units.c^2)+1).^2);
gamma = 1./sqrt(1-(v/units.c).^2);

nrows = 3;
ncols = 1;
isub = 1;

if 0
hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,E*1e-6,EK*1e-6)
hca.XLabel.String = 'Energy (MeV)';
hca.YLabel.String = 'Kinetic energy (MeV)';
end

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,E*1e-6,v*1e-3,E*1e-6,E*0+units.c*1e-3)
hca.XLabel.String = 'Energy (MeV)';
hca.YLabel.String = 'Speed (km/s)';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,E*1e-6,v/units.c)
hca.XLabel.String = 'Energy (MeV)';
hca.YLabel.String = 'Speed (c)';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,E*1e-6,gamma)
hca.XLabel.String = 'Energy (MeV)';
hca.YLabel.String = '\gamma';




%% Dipole field
units = irf_units;
RE = units.RE;
k0 = -8.02e15;

x = RE*linspace(-5,5,100);
y = 0;
z = RE*linspace(-5,5,100);

[X,Y,Z] = meshgrid(x,y,z);

[BX,BY,BZ,R,EL,AZ,POL] = magnetic_dipole_field(X,Y,Z,k0);

BX(R/RE<1) = NaN;
BY(R/RE<1) = NaN;
BZ(R/RE<1) = NaN;

nrows = 4;
ncols = 2;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,squeeze(X)/RE,squeeze(Z)/RE,squeeze(R)/RE)
hcb = colorbar('peer',hca); hcb.YLabel.String = 'R/R_E';
%axis(hca,'equal')

hca = h(isub); isub = isub + 1;
pcolor(hca,squeeze(X)/RE,squeeze(Z)/RE,squeeze(EL))
hcb = colorbar('peer',hca); hcb.YLabel.String = 'Elevation angle';

hca = h(isub); isub = isub + 1;
pcolor(hca,squeeze(X)/RE,squeeze(Z)/RE,squeeze(AZ))
hcb = colorbar('peer',hca); hcb.YLabel.String = 'Azimuthal angle';

hca = h(isub); isub = isub + 1;
pcolor(hca,squeeze(X)/RE,squeeze(Z)/RE,squeeze(POL))
hcb = colorbar('peer',hca); hcb.YLabel.String = 'Polar angle';


hca = h(isub); isub = isub + 1;
pcolor(hca,squeeze(X)/RE,squeeze(Z)/RE,squeeze(BX))
hcb = colorbar('peer',hca); hcb.YLabel.String = 'B_x';

hca = h(isub); isub = isub + 1;
pcolor(hca,squeeze(X)/RE,squeeze(Z)/RE,squeeze(BY))
hcb = colorbar('peer',hca); hcb.YLabel.String = 'B_y';

hca = h(isub); isub = isub + 1;
pcolor(hca,squeeze(X)/RE,squeeze(Z)/RE,squeeze(BZ))
hcb = colorbar('peer',hca); hcb.YLabel.String = 'B_z';


for ip = 1:nrows*ncols
  shading(h(ip),'flat')
end
%% Integrate particle trajectories
units = irf_units;
RE = units.RE;
E_eV = 0.2e6; % eV
m = units.mp;
q = units.e;
v = units.c*sqrt(1-1./(E_eV.*units.e./(m*units.c^2)+1).^2);
v_unit = [1, 0, 2]; v_unit = v_unit/norm(v_unit);
v_init = v*v_unit;
x_init = [-4, 0, 0]*RE;

options = odeset('RelTol',1e-9,'AbsTol',1e-9);

        
EoM = @(t,x) eom_dipole(t,x,m,q); 
[t,rv] = ode45(EoM,[0 60*10],[x_init,v_init],options);

%%
nrows = 4;
ncols = 1;
isub = 1;



if 1 % 3D
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  scatter3(rv(:,1)/RE,rv(:,2)/RE,rv(:,3)/RE,1,t/60)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'Time (min)';

  hold(hca,'on')
  plot3(hca,rv(1,1)/RE,rv(1,2)/RE,rv(1,3)/RE,'go',...
            rv(end,1)/RE,rv(end,2)/RE,rv(end,3)/RE,'rx')
  hold(hca,'off')
  %hcb = colorbar('peer',hca); 
  hca.ZLabel.String = 'x/R_E';
  hca.YLabel.String = 'y/R_E';
  hca.ZLabel.String = 'z/R_E';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on'; 
  hca.XLim = 2*[-10,10];
  hca.ZLim = 2*[-10,10];
  hca.YLim = 2*[-10,10];

  hold(hca,'on')
  plot3(hca,hca.XLim(2)+rv(:,1)*0,rv(:,2)/RE,rv(:,3)/RE   ,'color',[0.7 0.7 0.7])
  plot3(hca,rv(:,1)/RE,hca.YLim(2)+rv(:,1)*0,rv(:,3)/RE,'color',[0.7 0.7 0.7])
  plot3(hca,rv(:,1)/RE,rv(:,2)/RE,hca.ZLim(1)+rv(:,3)*0,'color',[0.7 0.7 0.7])
  hold(hca,'off')
end
if 1 % 2D (x,y)
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  scatter(rv(:,1)/RE,rv(:,2)/RE,1,t/60)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'Time (min)';
  hca.XLim = 2*[-10,10];
  hca.YLim = 2*[-10,10];
end
if 1 % (t,W)
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  plot(hca,t/60,units.me*sqrt(rv(:,4).^2 + rv(:,5).^2 + rv(:,6).^2)/2)
  hcb = colorbar('peer',hca);
  hcb.XLabel.String = 'Time (min)';  
  hcb.YLabel.String = 'Kinetic';  
end
if 1 % (t,z)
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  plot(hca,t/60,rv(:,3)/RE)
  hcb = colorbar('peer',hca);
  hcb.XLabel.String = 'Time (min)';
  hcb.YLabel.String = 'z (r_e)';
  
end





