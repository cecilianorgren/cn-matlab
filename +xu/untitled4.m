iOMNI = iPDist1.omni;
iEDensity = iOMNI;
iEDensity.data = iEDensity.data.*iEDensity.depend{1}.^2/units.mp/units.mp*4*units.e;
c_eval('WB? = gseB?.abs2*1e-18/2/units.e;',1:4)

hca = subplot(3,1,1);
irf_plot(hca,iEDensity.specrec('energy'))
hca.YScale = 'log';

hca = subplot(3,1,2);
irf_plot(hca,OMNI.deflux.specrec('energy'))
hca.YScale = 'log';

hca = subplot(3,1,3);
irf_plot(hca,WB1)

