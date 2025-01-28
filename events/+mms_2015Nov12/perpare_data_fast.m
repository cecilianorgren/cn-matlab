ic = 1:4;
% Current: curlometer
c_eval('gseR?fast = gseR?.resample(gseB?srvy);',1:4)
[Jcurlsrvy,divBbrstsrvy,Bbrstsrvy,JxBbrstsrvy,divTshearbrstsrvy,divPbbrstsrvy] = c_4_j('gseR?fast','gseB?srvy');
gseJcurlsrvy = irf.ts_vec_xyz(Jcurlsrvy.time,Jcurlsrvy.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurlsrvy.data = gseJcurlsrvy.data*1e9; Jcurlsrvy.units = 'nAm^{-2}';
gseJcurlsrvy.time = EpochTT(gseJcurlsrvy.time); gseJcurlsrvy.name = '4sc current density';

%% Calculate currents from moments
c_eval('gseJe?fast = -units.e*ne?fast*gseVe?fast*1e3*1e6*1e9; gseJe?fast.units = ''nA/m^2''; gseJe?fast.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi?fast = units.e*ne?fast*gseVi?fast.resample(ne?fast.time)*1e3*1e6*1e9; gseJi?fast.units = ''nA/m^2''; gseJi?fast.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ?fast = (gseJe?fast+gseJi?fast);',ic);
gseAvJfast = (gseJ1fast+gseJ2fast.resample(gseJ1fast.time)+gseJ3fast.resample(gseJ1fast.time)+gseJ4fast.resample(gseJ1fast.time))/4; 
