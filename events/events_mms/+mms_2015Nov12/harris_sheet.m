units = irf_units;
B0 = 10;
L = 1.2;
Te = 50;

Bx_Harris = @(z) -B0*tanh(z/L);
Bx_HarrisL = @(z,L) -B0*tanh(z/L);
Jy_Harris = @(z) -B0*1e-9*(1 - tanh(z/L).^2)/units.mu0/(L*1e3)*1e9;

tint_CS = irf.tint('2015-11-12T07:19:20.45Z/2015-11-12T07:19:21.80Z');
B_Obs = mvaAvB.tlim(tint_CS);
B_Offset = -2.3;
velocity = 70;
z_Obs = (B_Obs.time-B_Obs.time.start +- 0.5*(B_Obs.time.stop-B_Obs.time.start))*velocity;

J_Obs = mvaAvJ.tlim(tint_CS);
z_Obs_J = (J_Obs.time-J_Obs.time.start +- 0.5*(J_Obs.time.stop-J_Obs.time.start))*velocity;

% Threshold current based on electron thermal velocities
vth = cn_eV2v(Te,'eV');
Jvth = -2*5e6*units.e*vth*1e3*1e9; % nA/m^2

z = linspace(-L*3,L*3,100);

fontsize = 16;
hca = subplot(2,1,1);
plot(hca,z,Bx_Harris(z),z_Obs,B_Obs.x.data-B_Offset,z,Bx_HarrisL(z,1.2))
hca.YLim = [-15 15];
hca.YLabel.String = 'B_L(n) [nT]';
hca.XLabel.String = 'n [km]';
colors = hca.ColorOrder;
hca.ColorOrder = colors(1,:);
irf_legend(hca,{'Harris current sheet: B_L(n) = B_0 tanh(x/L)'},[0.98 0.98],'fontsize',fontsize)
irf_legend(hca,{['B_0 = ' num2str(B0) ', L = ' num2str(L) ]},[0.98 0.85],'fontsize',fontsize)
hca.ColorOrder = colors(2,:);
irf_legend(hca,{'Observed B_L'},[0.98 0.70],'fontsize',fontsize)
hca.FontSize = fontsize;

hca = subplot(2,1,2);
plot(hca,z,Jy_Harris(z),z_Obs_J,J_Obs.y.data,z,Jy_Harris(z)*0+Jvth)
