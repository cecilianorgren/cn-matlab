% Harris current sheet
tcenter = irf_time('2017-06-17T20:24:07.15','utc>epochtt');
tcenter = irf_time('2015-11-12T07:19:21.2','utc>epochtt');
tintObs = tcenter + [-1 1];
ic = 1;

c_eval('mvaB? = gseB?.resample(gseE1)*lmn''; mvaB?.name = ''B LMN'';',ic)
c_eval('mvaJ? = gseJ?*lmn''; mvaJ?.name = ''J LMN'';',ic)

c_eval(['obsB = mvaB?.tlim(tintObs);' ...
        'obsJ = mvaJ?.tlim(tintObs);'],ic)
zObsB = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsJ = (obsJ.time.epochUnix-mean(obsJ.time.epochUnix))*CS_normal_velocity;

CS_normal_velocity = 67;
L = 16;
dL = diff(zObsB(1:2));

hca = subplot(2,1,1);
plot(hca,zObsB,obsB.data,zObsB,-1-11*tanh(zObsB/L))
hca.YLabel.String = 'B (nT)';
hca.XLabel.String = 'N (km)';
legend(hca,'L','M','N','Fit L','location','best')
hca.Title.String = sprintf('MMS%g, half width = %g km',ic,L);

hca = subplot(2,1,2);
dB = diff(1+12*tanh(zObsB/L));
J = (dB*1e-9)/(dL*1e3)/units.mu0;
plot(hca,zObsJ,obsJ.data,zObsB(2:end)-0.5*dL,J*1e9)
hca.YLabel.String = 'J (nA/m^2)';
hca.XLabel.String = 'N (km)';
legend(hca,'L','M','N','Fit \Delta B_L/\Delta N','location','best')