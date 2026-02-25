% environment, density [cc], temperature [eV], magnetic field [nT]
environments = {'photosphere',1,1,1;...
                'chromosphere',1e10,1e0,1;... % mean free path, lam_c = 1e3 km
                'solar corona',1e7,1e3,1;... % mean free path, lam_c = 1 km
                'solar wind',1e1,1e1,1;... % mean free path, lam_c = 1e7 km
                'metals',1,1,1;...
                'ionosphere',1e6,1e-1,1;...
                'magnetosphere',1e1,1e3,1;...
                'fusion',1,1,1;...
                'center of Sun',1,1,1;...
                'radiation belts',1,1,1;...
                'interstellar medium',1,1,1;...
                'interplanetary',1,1,1;...
                'intergalactic',1,1,1,...
                };

    
n = [environments{:,2}];
T = [environments{:,3}];
B = [environments{:,4}];

hca = subplot(1,1,1);
loglog(hca,T,n,'o')
hca.XLabel.String = 'Temperature (eV)';
hca.YLabel.String = 'Density (cm^{-3})';

for ie = 1:size(environments,1)
  ht = text(T(ie),n(ie),environments{ie,1});
  ht.HorizontalAlignment = 'center';
end

hca.YLim = hca.YLim.*[1e-1 1e1];
hca.XLim = hca.XLim.*[1e-1 1e1];
  