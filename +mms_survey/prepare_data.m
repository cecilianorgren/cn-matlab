units = irf_units;

% Calculate curlometer current
disp('Calculating ExB velocity, speeds, currents, frequencies, length scales...')
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';

% Calculate currents from moments
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);
gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; 

% Calculate gradients
gseGradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4); gseGradPe.units = 'nPa/km';
gseGradPi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePi1,gsePi2,gsePi3,gsePi4); gseGradPi.units = 'nPa/km';
gseGradTe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gseTe1,gseTe2,gseTe3,gseTe4); gseGradTe.units = 'eV/km';
gseGradTi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gseTi1,gseTi2,gseTi3,gseTi4); gseGradTi.units = 'eV/km';
gseGradNe = c_4_grad('gseR?','ne?','grad');
gseGradNi = c_4_grad('gseR?','ni?','grad');

% Calculate convective electric fields
c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.units = ''mV/m'';',ic)
c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.units = ''mV/m'';',ic)

% Pitchangle distribution
%c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,18);',ic)

%% Wave analysis
if 0
  %%
c_eval('wavE? = irf_wavelet(gseE?.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
c_eval('wavE?.f_units = ''Hz''; wavE?.f_label = ''f [Hz]''; wavE?.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
c_eval('wavB? = irf_wavelet(gseB?scm.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
c_eval('wavB?.f_units = ''nT''; wavB?.f_label = ''f [Hz]''; wavB?.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)
%%
tintPol = irf.tint('2015-10-16T10:33:40.00Z/2015-10-16T10:33:52.00Z');
c_eval('polarization? = irf_ebsp(gseE?.tlim(tintPol),gseB?scm.tlim(tintPol),gseB?.tlim(tintPol),gseB?.tlim(tintPol),gseR?brsttime.tlim(tintPol),[10 2200],''polarization'',''fac'');',1)
end
%% Calculate some additional parameters
% ExB velocity
units = irf_units;
c_eval('gseVExB? = cross(gseE?,gseB?.resample(gseE?.time))/gseB?.abs.resample(gseE?.time)/gseB?.abs.resample(gseE?.time)*1e3;',ic) % km/s

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?);',ic)
 
% Speeds
c_eval('matB? = gseB?.abs.data,ne?.resample(gseB?.time).data;',ic)
c_eval('matParTe? = facTe?.xx.resample(gseB?.time).data;',ic)
c_eval('matParTi? = facTi?.xx.resample(gseB?.time).data;',ic)
c_eval('matPerTe? = (facTe?.yy.resample(gseB?.time).data + facTe?.zz.resample(gseB?.time).data)/2;',ic)
c_eval('matPerTi? = (facTi?.yy.resample(gseB?.time).data + facTi?.zz.resample(gseB?.time).data)/2;',ic)
c_eval('matTe? = facTe?.trace.resample(gseB?.time).data/3;',ic)
c_eval('matTi? = facTi?.trace.resample(gseB?.time).data/3;',ic)
c_eval('matNe? = ne?.resample(gseB?.time).data;',ic)

c_eval('vte?perp = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matTi?,''Vte''); vte?perp = irf.ts_scalar(gseB?.time,vte?perp)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vte?par = irf_plasma_calc(matB?,matNe?,0,matParTe?,matTi?,''Vte''); vte?par = irf.ts_scalar(gseB?.time,vte?par)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vte? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Vte''); vte? = irf.ts_scalar(gseB?.time,vte?)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vtp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matPerTi?,''Vtp''); vtp? = irf.ts_scalar(gseB?.time,vtp?)*1e-3; vtp?.units = ''km/s'';',ic)
c_eval('vA? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Va''); vA? = irf.ts_scalar(gseB?.time,vA?)*1e-3; vA?.units = ''km/s'';',ic)

% Frequencies
c_eval('flh? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Flh''); flh? = irf.ts_scalar(gseB?.time,flh?);',ic)
c_eval('fce? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fce''); fce? = irf.ts_scalar(gseB?.time,fce?);',ic)
c_eval('fcp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fcp''); fcp? = irf.ts_scalar(gseB?.time,fcp?);',ic)
c_eval('fpe? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fpe''); fpe? = irf.ts_scalar(gseB?.time,fpe?);',ic)
c_eval('fpp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fpp''); fpp? = irf.ts_scalar(gseB?.time,fpp?);',ic)

% Length scales
c_eval('Lp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Li''); Lp? = irf.ts_scalar(gseB?.time,Lp?)*1e-3; Lp?.units = ''km'';',ic)
c_eval('Le? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Le''); Le? = irf.ts_scalar(gseB?.time,Le?)*1e-3; Le?.units = ''km'';',ic)
c_eval('Ld? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Ld''); Ld? = irf.ts_scalar(gseB?.time,Ld?)*1e-3; Ld?.units = ''km'';',ic)
c_eval('re? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Roe''); re? = irf.ts_scalar(gseB?.time,re?)*1e-3; re?.units = ''km'';',ic)
c_eval('rp? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Rop''); rp? = irf.ts_scalar(gseB?.time,rp?)*1e-3; rp?.units = ''km'';',ic)

c_eval('beta? = (re?/Le?).^2;',ic)
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)

% JxB
c_eval('gseJxB? = gseJ?.cross(gseB?.resample(gseJ?));',ic)

% Electron pitchangles
c_eval('ePitch? = ePDist?.pitchangles(dmpaB?);',ic);

% Non-ideal electric field, E+VexB
c_eval('gseEVexB? = gseE?.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''E+VexB'';',ic)

c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('[gseVi?par,gseVi?perp] = irf_dec_parperp(gseB?,gseVi?); gseVi?par.name = ''Vi par''; gseVi?perp.name = ''Vi perp'';',ic)

% Parallel and perpendicular electric fields
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)


%% EDR signatures
c_eval('facPepp? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('facPeqq? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic);

% Compute Q and Dng from facPepp
c_eval('Q? = (facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2)./(facPepp?.yy.data.^2+2*facPepp?.yy.data.*facPepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
c_eval('Dng? = sqrt(8*(facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2))./(facPepp?.xx.data+2*facPepp?.yy.data);',ic);
c_eval('Dng? = irf.ts_scalar(ne?.time,Dng?);',ic);

% Compute agyrotropy Aphi from facPeqq
c_eval('agyro? = 2*(facPeqq?.yy-facPeqq?.zz)/(facPeqq?.yy+facPeqq?.zz); agyro? = agyro?.abs',ic);

% Compute temperature ratio An
c_eval('Temprat? = facPepp?.xx/(facPepp?.yy);',ic);

% Compute electron Mach number
c_eval('Me?perp = gseVe?perp.abs/vte?perp;',ic);
c_eval('Me?par = gseVe?par.abs/vte?par;',ic);

% Compute current density and J.E
c_eval('EdotJ? = gseE?.resample(gseJ?).dot(gseJ?)/1000; EdotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)

% Calculate epsilon and delta parameters
c_eval('oce? = fce?*2*pi;',ic)
c_eval('EdotVe? = gseE?.resample(gseVe?).dot(gseVe?);',ic);
c_eval('epsilone? = abs(6*pi*EdotVe?/(oce?.resample(gseVe?)*(facTe?.trace)));',ic);

c_eval('gseVexBE? = gseE?.resample(gseVe?)+gseVexB?.data;',ic);

c_eval('deltae? = gseVexB?/(gseVe?perp.abs*gseB?.resample(gseVe?).abs*1e-9);',ic);
%c_eval('deltae? = irf.ts_scalar(Uevec?.time,deltae?);',ic);

%% Average properties
avNe = (ne1+ne2.resample(ne1.time)+ne3.resample(ne1.time)+ne4.resample(ne1.time))/4; avNe.name = '<ne>';
gseAvE = (gseE1+gseE2.resample(gseE1.time)+gseE3.resample(gseE1.time)+gseE4.resample(gseE1.time))/4; 
gseAvVe = (gseVe1+gseVe2.resample(gseVe1.time)+gseVe3.resample(gseVe1.time)+gseVe4.resample(gseVe1.time))/4; 
gseAvVeperp = (gseVe1perp+gseVe2perp.resample(gseVe1perp.time)+gseVe3perp.resample(gseVe1perp.time)+gseVe4perp.resample(gseVe1perp.time))/4; 
gseAvVi = (gseVi1+gseVi2.resample(gseVi1.time)+gseVi3.resample(gseVi1.time)+gseVi4.resample(gseVi1.time))/4; 
gseAvB = (gseB1+gseB2.resample(gseB1.time)+gseB3.resample(gseB1.time)+gseB4.resample(gseB1.time))/4; 
gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; gseAvJ.units = gseJ1.units;
gseAvVExB = (gseVExB1+gseVExB2.resample(gseVExB1.time)+gseVExB3.resample(gseVExB1.time)+gseVExB4.resample(gseVExB1.time))/4; 
gseAvVexB = (gseVexB1+gseVexB2.resample(gseVexB1.time)+gseVexB3.resample(gseVexB1.time)+gseVexB4.resample(gseVexB1.time))/4; 
gseAvVixB = (gseVixB1+gseVixB2.resample(gseVixB1.time)+gseVixB3.resample(gseVixB1.time)+gseVixB4.resample(gseVixB1.time))/4; 
gseAvB = (gseB1+gseB2.resample(gseB1.time)+gseB3.resample(gseB1.time)+gseB4.resample(gseB1.time))/4;

% Ohm's law terms
%gseAvEVexB = (gseEVexB1 + gseEVexB2.resample(gseEVexB1) + gseEVexB3.resample(gseEVexB1) + gseEVexB4.resample(gseEVexB1))/4;
gseOhmGradPe = gseGradPe/avNe.resample(gseGradPe.time)/units.e*1e-9*1e-6; gseOhmGradPe.units = 'mV/m';
gseOhmVexB = gseAvVexB; gseOhmVexB.units = 'mV/m';
gseOhmVixB = gseAvVixB; gseOhmVixB.units = 'mV/m';
gseOhmJxB_a = gseAvJ.resample(gseAvB.time).cross(gseAvB)/avNe/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_a.units = 'mV/m';
gseOhmJxB_b = gseJcurl.resample(gseAvB.time).cross(gseAvB)/avNe/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_b.units = 'mV/m';
gseOhmJxB_c = (gseJxB1/ne1+gseJxB2.resample(gseJxB1.time)/ne2+gseJxB3.resample(gseJxB1.time)/ne3+gseJxB4.resample(gseJxB1.time)/ne4)/4/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB_c.units = 'mV/m'; 
gseOhmJxB = gseOhmJxB_c;

%% Misc
ic=1:4;

avPe = (gsePe1.trace+gsePe2.trace.resample(gsePe1)+gsePe3.trace.resample(gsePe1)+gsePe3.trace.resample(gsePe1))/3/4;
LgradP = avPe/gseGradPe.abs; LgradP.name = 'L_{P}';

% Proxy for magnetic field slippage (Scudder2105)
avFce = (fce1+fce2.resample(fce1)+fce3.resample(fce1)+fce4.resample(fce1))/4;
avRhoe = (re1+re2.resample(re1)+re3.resample(re1)+re4.resample(re1))/4;
Y = avFce*(avRhoe/LgradP).^2; Y.units = 's^-1';
Ynodim = (avRhoe/LgradP); Y.units = 's^-1';

facAvTe = (facTe1+facTe2.resample(facTe1)+facTe3.resample(facTe1)+facTe4.resample(facTe1))/4;

% Curl of E+VexB (Scudder2105)
gseRotRe = mms_2015Oct16.rotRe(gseR1,gseR2,gseR3,gseR4,gseEVexB1,gseEVexB2,gseEVexB3,gseEVexB4);
