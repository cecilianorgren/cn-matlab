%% 2d map of large drift vs small drift regime, vExV/vti
units = irf_units;
Ti = linspace(1000,7000,100); % eV
B = 10e-9;
E = linspace(0,30,100)*1e-3;

[Ti_,E_] = meshgrid(Ti,E);

vExB = E_/B*1e-3; % km/s
vTi = sqrt(Ti_*units.eV*2/units.mp)*1e-3; % km/s

nrows = 2;
ncols = 2;

isub = 1;
hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,Ti,E*1e3,vExB./vTi);
shading(hca,'flat')
colormap(hca,cn.cmap('bluered3'));
hcb = colorbar('peer',hca);
hca.YLabel.String = sprintf('E (mV/m), B = %g nT',B*1e9);
hca.XLabel.String = 'T_i (eV)';
hca.Title.String = 'v_{ExB}/v_{Ti}';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,vTi,E/B*1e-3,vExB./vTi);
shading(hca,'flat')
colormap(hca,cn.cmap('bluered3'));
hcb = colorbar('peer',hca);
hca.YLabel.String = sprintf('E/B (km/s)');
hca.XLabel.String = 'v_{Ti} (km/s)';
hca.Title.String = 'v_{ExB}/v_{Ti}';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,Ti,E*1e3,log10(vExB./vTi));
shading(hca,'flat')
colormap(hca,cn.cmap('bluered3'));
hcb = colorbar('peer',hca);
hca.CLim = [-1 1];
hca.YLabel.String = sprintf('E (mV/m), B = %g nT',B*1e9);
hca.XLabel.String = 'T_i (eV)';
hca.Title.String = 'log_{10}(v_{ExB}/v_{Ti})';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,vTi,E/B*1e-3,log10(vExB./vTi));
shading(hca,'flat')
colormap(hca,cn.cmap('bluered3'));
hca.CLim = [-1 1];
hcb = colorbar('peer',hca);
hca.YLabel.String = sprintf('E/B (km/s)');
hca.XLabel.String = 'v_{Ti} (km/s)';
hca.Title.String = 'log_{10}(v_{ExB}/v_{Ti})';

%% Density fluctuations in wave potential
Phi = logspace(1,4,10);
Phi = linspace(700,1200,20);
Te = logspace(2.5,3.3,20);
Te = linspace(500,2000,20);
n0 = 1;

[Te_,Phi_] = meshgrid(Te,Phi);
ne = n0*exp(Phi_./Te_);

nrows = 1;
ncols = 2;

isub = 1;
hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,Te,Phi,ne);
shading(hca,'flat')
colormap(hca,cn.cmap('bluered3'));
hcb = colorbar('peer',hca);
hca.XScale = 'lin';
hca.YScale = 'lin';
hca.YLabel.String = '\Phi (V)';
hca.XLabel.String = 'T_e (eV)';
hca.Title.String = 'n_e(\Phi) = n_0exp(e\Phi/k_BT_e)';
hca.CLim = [1 4];
%hca.XTick = 0:100:Te(end);

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,Te,Phi,log10(ne));
shading(hca,'flat')
colormap(hca,cn.cmap('bluered3'));
hcb = colorbar('peer',hca);
hca.XScale = 'lin';
hca.YScale = 'lin';
hca.XMinorTick = 'on';
hca.YLabel.String = '\Phi (V)';
hca.XLabel.String = 'T_e (eV)';
hca.Title.String = 'log_{10} n_e(\Phi) = log_{10} n_0exp(e\Phi/k_BT_e)';

%% Input distributions 
units = irf_units;
n = 0.05;
Ti = 4000; % eV
Te = 1000;
B = 10e-9;
E = 25e-3;


vExB = E/B*1e-3; % km/s
vTi = sqrt(Ti*units.eV*2/units.mp)*1e-3; % km/s
vTe = sqrt(Te*units.eV*2/units.me)*1e-3; % km/s
v = linspace(-vTi,0.5*vTe,1000);

ve = vExB*(1+Te/Ti); % second term is electron diamagnetic drift
vi = 0;
                    
fe = cn.maxwellian(v,Te,n,ve,'e');
fi = cn.maxwellian(v,Ti,n,vi,'p');

nrows = 1;
ncols = 1;

isub = 1;
hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,v,[fe, fi],'-',...         
         vExB*[1 1],[0 max(fi)],'-',...
         ve*[1 1],[0 max(fi)],'-',...
         vTi*[1 1],[0 max(fi)],'-')

legend('f_e','f_i','v_{ExB}','v_{de}','v_{Ti}')
hca.Title.String = sprintf('v_{ExB}/v_{Ti} = %g',vExB/vTi);
