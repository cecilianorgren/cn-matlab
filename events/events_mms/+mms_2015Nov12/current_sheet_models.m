units = irf_units;
B0 = 1;
L = 1;
Te = 50;

% Currnet sheet models
% Harris current sheet
Bx_Harris = @(z,L) B0*tanh(z/L);
Jy_Harris = @(z,L) B0*(1 - tanh(z/L).^2)/(L);
% Janaki2012 current sheet (taken from Yoon2014, because Janaki had typos)
Bx_Janaki = @(z,r) (1+0.5*r)^0.5*(1+r)*sinh(sqrt(1+0.5*r).*z).*cosh(sqrt(1+0.5*r).*z)./((1+r)*cosh(sqrt(1+0.5*r).*z).^2-0.5*r);
Jy_Janaki = @(x,r) (1+0.5*r)*(1+r)*(cosh(1*sqrt(1+0.5*r).*z).^2+0.5*r)./(((1+r).*cosh(sqrt(1+0.5*r).*z).^2-0.5*r).^2);

% Yoon2014 current sheet, Bx and Jy remains unchanged w.r.t. Janaki's
Bx_Janaki = @(z,r) (1+0.5*r)^0.5*(1+r)*sinh(sqrt(1+0.5*r).*z).*cosh(sqrt(1+0.5*r).*z)./((1+r)*cosh(sqrt(1+0.5*r).*z).^2-0.5*r);
Jy_Janaki = @(x,r) (1+0.5*r)*(1+r)*(cosh(1*sqrt(1+0.5*r).*z).^2+0.5*r)./(((1+r).*cosh(sqrt(1+0.5*r).*z).^2-0.5*r).^2);




r = -0.88;
z = linspace(-L*5,L*5,100);
nRows = 4;
hca = subplot(nRows,1,1);
plot(hca,z,Bx_Janaki(z,r),z,Bx_Harris(z,L))

hca = subplot(nRows,1,2);
plotyy(hca,z,Jy_Janaki(z,r),z,Jy_Harris(z,L))

hca = subplot(nRows,1,3);
plot(hca,z,cosh(z),z,sinh(z),z,tanh(z))

hca = subplot(nRows,1,4);
plotyy(hca,z,Jy_Janaki_top(z,r),z,Jy_Janaki_bot(z,L))

%%
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
